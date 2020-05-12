using Plots
include("heat_equation.jl")


# Heat Equation
# Case 1:
# Constant heat sources on each end of an insulated rod.
function heatequation1()
    rod = HeatEquation.HeatEquation1D(
        L = 10,
        leftBoundary = 1,
        rightBoundary = 1,
        T = 30
    )
    solution = HeatEquation.solve(rod)
    anim = HeatEquation.plot(rod, solution, 50)
    gif(anim, "./plots/heat_equation_case_1.gif")
end


# Case 2:
# Both ends of an insulated rod are held at absolute zero, with uniform initial condition.
function heatequation2()
    rod = HeatEquation.HeatEquation1D(
        L = 10,
        leftBoundary = 0,
        rightBoundary = 0,
        ic = x->1,
        T = 30
    )
    solution = HeatEquation.solve(rod)
    anim = HeatEquation.plot(rod, solution, 50)
    gif(anim, "./plots/heat_equation_case_2.gif")
end


# Case 3:
# Non-uniform Dirichlet boundary conditions, uniform initial conditions.
function heatequation3()
    rod = HeatEquation.HeatEquation1D(
        L = 10,
        leftBoundary = 0,
        rightBoundary = 1,
        T = 30
    )
    solution = HeatEquation.solve(rod)
    anim = HeatEquation.plot(rod, solution, 50)
    gif(anim, "./plots/heat_equation_case_3.gif")
end

heatequation1()
heatequation2()
heatequation3()