# Numerical Methods for Conservation Laws: Project 1

Welcome to the repository for the 1st project of EPFL MATH-459: Numerical Methods for Conservation Laws. This project is dedicated to the implementation of numerical simulations for a specific Partial Differential Equation (PDE) problem, that is the shallow water equations, as outlined in the assignment.

## Code Structure

The project code is organized into three main files:

1. `part1_1.m`
2. `part1_2.m`
3. `part1_3.m`

Each of these files executes a numerical simulation for the PDE problem specified in the project assignment. It is important to note that the user can control whether to visualize an animation of the solution by modifying the first line of each file, setting `animation` to either "true" or "false."

## Core Functions

### 1. `conservative_scheme.m`

The generic numerical conservative scheme is implemented in this function. It takes several inputs, including a physical flux and a numerical flux, and plays a central role in the numerical simulations conducted in the project.

### 2. `flux_phys.m`

The physical flux for the shallow water equations, a key component of the conservative scheme, is implemented in this file.

### 3. Numerical Flux Implementations

Two different numerical flux methods, namely Lax-Friedrichs and Lax-Wendroff, are implemented and available in their respective files with the following names:

- `Lax_Friedrichs_flux.m`
- `Lax_Wendroff_flux.m`

## Usage

To execute a specific part of the project, run the corresponding `part1_X.m` file. Adjust the `animation` variable in the first lines to control the visualization of the solution.

## Project Contributors

This project was realized by:

- Francesco Sala
- Nicolo' Viscusi

Completed in December 2023.

If you have any questions or suggestions, please don't hesitate to reach out to the contributors.

