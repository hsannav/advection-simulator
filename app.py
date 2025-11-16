import streamlit as st
import numpy as np
import plotly.graph_objects as go

st.set_page_config(layout="wide", page_title="Advection PDE Simulator")

st.title("One-Dimensional Time-Dependent Advection Equations")

with st.sidebar:
    st.header("Simulation Parameters")
    
    wave_shape = st.selectbox(
        "Initial Condition Shape",
        ["Gaussian Bell", "Square Pulse", "Triangle Wave", "Cosine Hat"]
    )

    L = st.slider("Domain Length ($L$)", 1.0, 100.0, 1.0)
    v = st.slider("Advection Velocity ($v$)", 0.1, 5.0, 1.0)
    T_max = st.slider("Simulation Time ($T$)", 0.1, 5.0, 1.0)
    
    Nx = st.select_slider("Spatial Resolution ($N_x$)", options=[50, 100, 200, 400, 800], value=100)
    
    C = st.slider("Courant Number ($C = v \\frac{dt}{dx}$)", 0.1, 2.0, 0.8)
    dx = L / float(Nx)
    dt = C * dx / v
    Nt = int(round(T_max / dt))
    st.write(f"$dt = {dt:.4f}$")
    st.write(f"$dx = {dx:.4f}$")
    st.write(f"$N_t = {Nt}$")
    st.markdown("---")

    st.subheader("Schemes to Compare")
    #show_exact = st.checkbox("Exact Solution", value=True)
    show_fecs = st.checkbox("Forward Euler Centered", value=True)
    show_upwind = st.checkbox("Upwind", value=True)
    show_leapfrog = st.checkbox("Leapfrog", value=True)
    show_lw = st.checkbox("Lax-Wendroff", value=True)

    st.subheader("Shape Parameters")
    sigma = st.slider("Sigma (Gaussian/Cosine)", 0.01, 0.5, 0.05)
    x0 = st.slider("Initial Position ($x_0$)", 0.0, L, L/4.0)

x = np.linspace(0, L, Nx + 1)
x = x[:-1]

def get_initial_condition(x, shape, sig, center, length):
    if shape == "Gaussian Bell":
        return np.exp(-0.5 * ((x - center) / sig)**2)
    elif shape == "Square Pulse":
        return np.where(np.abs(x - center) < sig, 1.0, 0.0)
    elif shape == "Triangle Wave":
        return np.where(np.abs(x - center) < sig, 1.0 - np.abs(x - center)/sig, 0.0)
    elif shape == "Cosine Hat":
        r = np.abs(x - center)
        cond = r < sig
        return np.where(cond, np.cos(np.pi * r / (2 * sig))**2, 0.0)
    return np.zeros_like(x)

u0 = get_initial_condition(x, wave_shape, sigma, x0, L)

solvers = {}
integrals = {}

if show_fecs:
    u = u0.copy()
    u_all = [u.copy()]
    mass = [np.sum(u) * dx]
    for n in range(Nt):
        u_next = u - 0.5 * C * (np.roll(u, -1) - np.roll(u, 1))
        u = u_next
        u_all.append(u.copy())
        mass.append(np.sum(u) * dx)
    solvers["FECS"] = u_all
    integrals["FECS"] = mass

if show_upwind:
    u = u0.copy()
    u_all = [u.copy()]
    mass = [np.sum(u) * dx]
    for n in range(Nt):
        u_next = u - C * (u - np.roll(u, 1))
        u = u_next
        u_all.append(u.copy())
        mass.append(np.sum(u) * dx)
    solvers["Upwind"] = u_all
    integrals["Upwind"] = mass

if show_leapfrog:
    u_prev = u0.copy()
    u_curr = u_prev - 0.5 * C * (np.roll(u_prev, -1) - np.roll(u_prev, 1))
    u_all = [u_prev.copy(), u_curr.copy()]
    mass = [np.sum(u_prev) * dx, np.sum(u_curr) * dx]
    
    for n in range(Nt - 1):
        u_next = u_prev - C * (np.roll(u_curr, -1) - np.roll(u_curr, 1))
        u_prev = u_curr
        u_curr = u_next
        u_all.append(u_curr.copy())
        mass.append(np.sum(u_curr) * dx)
    solvers["Leapfrog"] = u_all
    integrals["Leapfrog"] = mass

if show_lw:
    u = u0.copy()
    u_all = [u.copy()]
    mass = [np.sum(u) * dx]
    for n in range(Nt):
        term1 = 0.5 * C * (np.roll(u, -1) - np.roll(u, 1))
        term2 = 0.5 * C**2 * (np.roll(u, -1) - 2*u + np.roll(u, 1))
        u_next = u - term1 + term2
        u = u_next
        u_all.append(u.copy())
        mass.append(np.sum(u) * dx)
    solvers["Lax-Wendroff"] = u_all
    integrals["Lax-Wendroff"] = mass


fig_anim = go.Figure()

colors = {"Exact": "black", "FECS": "red", "Upwind": "blue", "Leapfrog": "green", "Lax-Wendroff": "orange"}
dash = {"Exact": "solid", "FECS": "dash", "Upwind": "dot", "Leapfrog": "dashdot", "Lax-Wendroff": "longdash"}
for name, data in solvers.items():
    fig_anim.add_trace(go.Scatter(
        x=x, 
        y=data[0], 
        mode="lines", 
        name=name,
        line=dict(color=colors.get(name, "gray"), dash=dash.get(name, "solid"))
    ))
frames = []
step_skip = max(1, Nt // 50)

for k in range(0, Nt + 1, step_skip):
    frame_data = []
    for name, data in solvers.items():
        frame_data.append(go.Scatter(x=x, y=data[k]))
    frames.append(go.Frame(data=frame_data, name=str(k)))
fig_anim.update(frames=frames)

fig_anim.update_layout(
    xaxis=dict(range=[0, L], title="Position (x)"),
    yaxis=dict(range=[-2, 2], title="u(x,t)"),
    updatemenus=[dict(
        type="buttons",
        buttons=[dict(label="Play", method="animate", args=[None, dict(frame=dict(duration=50, redraw=True), fromcurrent=True)])]
    )]
)

st.plotly_chart(fig_anim, use_container_width=True)

if len(solvers) > 0:
    for solver in solvers.keys():
        data_matrix = np.array(solvers[solver])
        t_vals_heat = np.linspace(0, T_max, Nt + 1)
        
        fig_heatmap = go.Figure(data=go.Heatmap(
            z=data_matrix,
            x=x,
            y=t_vals_heat,
            colorscale='Viridis',
            colorbar=dict(title="u(x,t)")
        ))

        fig_heatmap.update_layout(
            title=f"Space-Time Diagram: {solver}",
            xaxis_title="Position (x)",
            yaxis_title="Time (t)",
            height=400
        )

        st.plotly_chart(fig_heatmap, use_container_width=True)