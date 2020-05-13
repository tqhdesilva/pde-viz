module WaveEquation
using DifferentialEquations
using DiffEqOperators
using LinearAlgebra

Base.@kwdef struct WaveEquation1D
    icₓ::Function
    ic::Function
    leftBoundary::Float64
    rightBoundary::Float64
    c²::Float64 = 1
    T::Float64
    L::Float64 = 10
    Δx::Float64 = 0.1
    Δt::Float64 = 0.005
end

function solve(d::WaveEquation1D)
    x = d.Δx:d.Δx:d.L - d.Δx
    t = d.Δt:d.Δt:d.T - d.Δt
    Nₓ = length(x)
    Nₜ = length(t)
    u₀ = d.ic.(x)
    du₀ = d.icₓ.(x)
    U₀ = [u₀; du₀]
    A₁₁ = zeros(Nₓ, Nₓ)
    A₁₂ = Matrix{Float64}(I, Nₓ, Nₓ)
    A₂₂ = zeros(Nₓ, Nₓ)

    A₂₁ = zeros(Nₓ, Nₓ)
    A₂₁[1, :] = [[-2, 1]; [0 for i = 1:Nₓ - 2]]
    A₂₁[end, :] = [[0 for i = 1:Nₓ - 2]; [1, -2]]
    for i = 2:Nₓ - 1
        A₂₁[i, :] = [
            [0 for j = 1:i - 2];
            [1, -2, 1];
            [0 for j = 1:Nₓ - i - 1]
        ]
    end
    A₂₁ = A₂₁ .* (1 / d.Δx^2) .* d.c²

    A = [
        [A₁₁ A₁₂];
        [A₂₁ A₂₂]
    ]

    function waveequation!(du, u, p, t)
        du[:] = A * u + [
            [d.leftBoundary];
            [0 for i = 1:Nₓ - 2];
            [d.rightBoundary];
            [0 for i = 1:Nₓ]
        ]
    end

    prob = ODEProblem(waveequation!, U₀, (0.0, d.T))
    sol = DifferentialEquations.solve(prob, KenCarp4())
    solarray = zeros(Nₓ, Nₜ + 1)
    solarray[:, 1] = u₀
    for (i, tᵢ) in enumerate(t)
        Uᵢ = sol(tᵢ)
        uᵢ = Uᵢ[1:Integer(length(Uᵢ) / 2)]
        solarray[:, i + 1] = uᵢ
    end
    return solarray
end
end