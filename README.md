# One-Dimensional Time-Dependent Advection Equations Simulator

## Description

This project is an interactive web application built with Python and Streamlit that simulates the one-dimensional time-dependent advection equation. It serves as a pedagogical tool to visualize and compare various finite-difference numerical schemes against the exact analytical solution.

The application allows users to modify simulation parameters in real-time—such as the Courant number, spatial resolution, and wave velocity—to observe how different discretization methods handle stability, dissipation, and dispersion.

## Mathematical Background

The simulator solves the linear advection equation, a partial differential equation (PDE) that describes the transport of a quantity (scalar field $u$) by a flow with velocity $v$:

$$
\frac{\partial u}{\partial t} + v \frac{\partial u}{\partial x} = 0
$$

Where:
* $u(x,t)$ is the scalar field being advected.
* $v$ is the constant velocity.

## Features

* **Interactive Simulation Controls:** Adjust Domain Length ($L$), Advection Velocity ($v$), Simulation Time ($T$), and Spatial Resolution ($N_x$).
* **Courant Number Analysis:** Dynamically adjust the Courant-Friedrichs-Lewy (CFL) condition ($C = v \frac{dt}{dx}$) to observe stability limits.
* **Multiple Initial Conditions:**
    * Gaussian Bell
    * Square Pulse
    * Triangle Wave
    * Cosine Hat
* **Numerical Schemes:** The application implements and compares the following schemes:
    * Exact Solution (Analytical)
    * Forward Euler Centered Space (FECS)
    * Upwind Scheme
    * Leapfrog Scheme
    * Lax-Wendroff Scheme
    * Crank-Nicolson (Implicit)
    * Backward Euler (Implicit)
* **Visualizations:**
    * Animated wave propagation (Plotly).
    * Space-Time (Hovmöller) heatmaps for every solver.
    * Mass conservation integrals to analyze numerical dissipation and instability.

## Installation and Usage

### Prerequisites

Ensure you have Python 3.8+ installed.

### Setup

1. Clone the repository:
   ```bash
   git clone https://github.com/hsannav/advection-simulator.git
   cd advection-simulator
   ```

2. Create a virtual environment (recommended):
   ```bash
   python -m venv venv
   source venv/bin/activate  # On Windows use `venv\Scripts\activate`
   ```

3. Install the required dependencies:
   ```bash
   pip install -r requirements.txt
   ```

4. Run the Streamlit application:
   ```bash
   streamlit run app.py
   ```

## Project Structure

* `app.py`: The main application script containing the PDE solvers and Streamlit interface.
* `requirements.txt`: List of Python dependencies.
* `README.md`: Project documentation.

## Numerical Methods Overview

This project demonstrates the behavior of different discretization strategies:

* **FECS:** Unconditionally unstable.
* **Upwind:** First-order accurate, stable if $C \leq 1$, but diffusive.
* **Leapfrog:** Second-order accurate, non-dissipative, but dispersive.
* **Lax-Wendroff:** Second-order accurate, minimizes dissipation compared to Upwind.
* **Implicit Schemes (CN, BE):** Unconditionally stable, requiring matrix inversion at each time step.

## Authors

* **Fernando Blanco**
* **Hugo Sánchez**

## License

This project is open-source and available for educational purposes.
