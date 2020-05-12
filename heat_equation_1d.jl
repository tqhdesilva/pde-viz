using Plots

Base.@kwdef struct HeatEquation1D
    Δx::Float64 = 0.1
    L::Float64
    leftBoundary::Float64
    rightBoundary::Float64
    ic::Function = x->0
    T::Float64 = 10
    Δt::Float64 = 0.005
    c2::Float64 = 1
end

function solve(d::HeatEquation1D)::Array{Float64,2}
    x = d.Δx:d.Δx:d.L - d.Δx
    t = d.Δt:d.Δt:d.T - d.Δt
    Nx = length(x)
    Nt = length(t)
    u0 = d.ic.(x)
    un = copy(u0)
    un[1] = d.leftBoundary
    un[end] = d.rightBoundary
    u = zeros(Nx, Nt + 1)
    u[:, 1] = u0
    u[:, 2] = un

    h = d.Δx / d.c2
    r = d.Δt / h^2
    A = zeros(Nx, Nx)
    A[1, 1] = 1
    A[end, end] = 1
    for i = 2:Nx - 1
        A[i, :] = [
            [0 for j = 1:i - 2];
            [r, 1 - 2 * r, r];
            [0 for j = 1:Nx - i - 1]
        ]
    end

    for i in 3:Nt + 1
        un = A * un
        u[:, i] = un
    end
    return u
end


function plot(d::HeatEquation1D, solution::Array{Float64,2}, every::Integer = 1)
    x = d.Δx:d.Δx:d.L - d.Δx
    t = d.Δt:d.Δt:d.T - d.Δt
    Nt = length(t)
    @gif for i = 1:Nt
        Plots.plot(x, solution[:, i], ylims = (minimum(solution), maximum(solution)))
    end when i % every == 0
end
