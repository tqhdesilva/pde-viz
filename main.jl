using Plots
include("heat_equation.jl")
include("wave_equation.jl")

function plot(
    d::Union{HeatEquation.HeatEquation1D,WaveEquation.WaveEquation1D},
    solution::Array{Float64,2},
    every::Integer = 1
)::Plots.Animation
    x = d.Δx:d.Δx:d.L - d.Δx
    t = d.Δt:d.Δt:d.T - d.Δt
    Nt = length(t)
    anim = @animate for i = 1:Nt
        Plots.plot(x, solution[:, i], ylims = (minimum(solution), maximum(solution)))
    end when i % every == 0
    return anim
end

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
    anim = plot(rod, solution, 50)
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
    anim = plot(rod, solution, 50)
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
    anim = plot(rod, solution, 50)
    gif(anim, "./plots/heat_equation_case_3.gif")
end


function waveequation1()
    wave = WaveEquation.WaveEquation1D(
        icₓ = x->0.0,
        ic = x->sin(2.0 * π * x / 10.0),
        leftBoundary = 0.0,
        rightBoundary = 0.0,
        L = 10.0,
        T = 10.0,
        Δt = 0.01
    )
    sol = WaveEquation.solve(wave)
    anim = plot(wave, sol, 10)
    gif(anim, "./plots/wave_equation_case_1.gif")
end

# TODO Convert to a first order system
# and apply stiff solver.
# This seems to be a stiff ODE.
# There's a bug in DifferentialEquations.jl
# that makes solving Dynamical ODEs(e.g. 2nd order ODE)
# impossible.
function waveequation2()
    wave = WaveEquation.WaveEquation1D(
        icₓ = x->-x,
        ic = x->x,
        leftBoundary = 0.0,
        rightBoundary = 0.0,
        L = 1.0,
        Δx = 0.005,
        T = 1.0,
        Δt = 0.00005,
        c² = 4
    )
    sol = WaveEquation.solve(wave)
    anim = plot(wave, sol, 200)
    gif(anim, "./plots/wave_equation_case_2.gif")
end

heatequation1()
heatequation2()
heatequation3()
waveequation1()
waveequation2()