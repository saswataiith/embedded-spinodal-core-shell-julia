# Fixed-particle embedded Cahn–Hilliard model, recovered July 2022 source.
# All quantities are dimensionless. The original CUDA sources are unchanged.
using FFTW, Random, Statistics, DelimitedFiles
Base.@kwdef struct Parameters
    A::Float64 = 2.0
    B::Float64 = 4.0
    chi::Float64 = 0.5
    P::Float64 = 5.0
    ca::Float64 = 0.5
    cb::Float64 = 0.0
    cc::Float64 = 1.0
    kappa::Float64 = 2.0
    mobility::Float64 = 1.0
    Q::Float64 = 0.5
end

function make_grid(n, spacing)
    k = collect(2pi .* FFTW.fftfreq(n, 1 / spacing))
    derivative_k = copy(k)
    iseven(n) && (derivative_k[n ÷ 2 + 1] = 0.0)
    kx = reshape(derivative_k, n, 1)
    ky = reshape(derivative_k, 1, n)
    k2 = reshape(k.^2, n, 1) .+ reshape(k.^2, 1, n)
    return (; kx, ky, k2, spacing)
end

function initial_fields(n; spacing=1.0, radius=70.0, width=4.0, seed=5749)
    # Radius and width are physical nondimensional lengths, not grid indices.
    phi = [0.5 * (1 - tanh((hypot((i-1-n÷2)*spacing,
            (j-1-n÷2)*spacing)-radius)/width)) for i in 1:n, j in 1:n]
    inside = phi .>= 0.5
    noise = 0.01 .* (2 .* rand(MersenneTwister(seed), n, n) .- 1)
    noise .-= mean(noise[inside])
    c = fill(0.5, n, n)
    c[inside] .+= noise[inside]
    return c, phi
end

function geometry(phi, p)
    h = @. phi^3 * (10 - 15phi + 6phi^2)
    g = @. phi^2 * (1-phi)^2
    return (; h, g, M=p.mobility .* h)
end

function bulk_energy(c, shape, p)
    h, g = shape.h, shape.g
    return @. p.A*(1-h)*(c-p.ca)^2 +
        p.B*h*(c-p.cb)^2*(c-p.cc)^2 + p.P*(1-p.chi*c)*g
end

function bulk_derivative(c, shape, p)
    h, g = shape.h, shape.g
    return @. 2p.A*(1-h)*(c-p.ca) +
        2p.B*h*(c-p.cb)*(c-p.cc)*(2c-p.cb-p.cc) - p.chi*p.P*g
end

function chemical_potential(c, shape, grid, p)
    return bulk_derivative(c, shape, p) .+
        real.(ifft(2p.kappa .* grid.k2 .* fft(c)))
end

function transport_rate(c, shape, grid, p)
    mu_hat = fft(chemical_potential(c, shape, grid, p))
    grad_x = real.(ifft(im .* grid.kx .* mu_hat))
    grad_y = real.(ifft(im .* grid.ky .* mu_hat))
    divergence_hat = im .* grid.kx .* fft(shape.M .* grad_x) .+
                     im .* grid.ky .* fft(shape.M .* grad_y)
    dissipation = sum(shape.M .* (grad_x.^2 .+ grad_y.^2)) * grid.spacing^2
    return divergence_hat, dissipation
end

function step(c, shape, grid, p, dt)
    rate, _ = transport_rate(c, shape, grid, p)
    denominator = @. 1 + 2p.Q*p.kappa*grid.k2^2*dt
    return real.(ifft(fft(c) .+ dt .* rate ./ denominator))
end

function energy(c, shape, grid, p)
    # Parseval form includes the even-grid Nyquist contribution to |grad c|².
    gradient_energy = p.kappa * sum(grid.k2 .* abs2.(fft(c))) / length(c)
    return (sum(bulk_energy(c, shape, p)) + gradient_energy) * grid.spacing^2
end

function run_simulation(c0, phi; dt=0.01, steps=1000, spacing=1.0,
                        p=Parameters(), directory=nothing)
    grid = make_grid(size(c0,1), spacing)
    shape = geometry(phi, p)
    c = copy(c0)
    history = zeros(steps+1, 6)
    for s in 0:steps
        _, dissipation = transport_rate(c, shape, grid, p)
        history[s+1,:] = [s*dt, mean(c), energy(c,shape,grid,p), minimum(c), maximum(c), dissipation]
        all(isfinite, c) || error("Non-finite composition at step $s")
        s < steps && (c = step(c, shape, grid, p, dt))
    end
    if directory !== nothing
        mkpath(directory)
        writedlm(joinpath(directory,"history.csv"), history, ',')
        writedlm(joinpath(directory,"composition.csv"), c, ',')
        writedlm(joinpath(directory,"particle.csv"), phi, ',')
    end
    return c, history
end

if abspath(PROGRAM_FILE) == @__FILE__
    c, phi = initial_fields(256)
    run_simulation(c, phi; directory=joinpath(@__DIR__, "results-256"))
    println("Saved 256 × 256 fields and diagnostics; final time = 10.")
end
