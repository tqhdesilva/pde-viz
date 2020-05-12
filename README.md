# Visualization for Partial Differential Equations

I made this repository to practice using Julia for solving and visualizing
partial differential equations. Right now there is only code for solving the
heat equation in 1D with Dirichlet boundary conditions.

After getting your Julia environment set up, run
```bash
julia --project=. main.jl
```
to reproduce the plots.


## Heat Equation
### Case 1
Uniform 0 initial conditions, with both ends of the rod held at 1.
![Case 1](./plots/heat_equation_case_1.gif)

### Case 2
The initial condition is 1 everywhere, with the boundaries set to 0.
![Case 2](./plots/heat_equation_case_2.gif)

### Case 3
0 initial condition, with one boundary condition at zero and the other at 1.
![Case 3](./plots/heat_equation_case_3.gif)

## Wave Equation
### Case 1
Standing sine wave.

![Case 1](./plots/wave_equation_case_1.gif)