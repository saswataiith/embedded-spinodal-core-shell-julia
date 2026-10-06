include("core_shell.jl")
using Test
p = Parameters()
c, phi = initial_fields(64; radius=16.0)
grid = make_grid(64, 1.0)
shape = geometry(phi,p)
v = randn(MersenneTwister(8), size(c))
eps = 1e-6
numerical = (energy(c .+ eps .* v,shape,grid,p)-energy(c .- eps .* v,shape,grid,p))/(2eps)
analytical = sum(chemical_potential(c,shape,grid,p) .* v)
@test isapprox(numerical, analytical; rtol=1e-6, atol=1e-6)
rate, dissipation = transport_rate(c,shape,grid,p)
@test isapprox(sum(chemical_potential(c,shape,grid,p) .* real.(ifft(rate))), -dissipation; rtol=1e-10)
@test abs(sum(real.(ifft(rate)))) < 1e-10
println("Energy derivative and discrete dissipation identity passed.")
finals = []
for dt in (0.01, 0.005, 0.0025)
    final, hist = run_simulation(c,phi; dt, steps=round(Int,10/dt),directory=joinpath(@__DIR__,"dt-$dt"))
    @test maximum(abs.(hist[:,2] .- mean(c))) < 1e-12
    @test maximum(diff(hist[:,3])) < 1e-8
    push!(finals,final)
    println("dt=$dt mass drift=$(maximum(abs.(hist[:,2].-mean(c)))) energy=$(hist[1,3]) → $(hist[end,3]) range=$(extrema(final))")
end
println("Final-field RMS dt .01/.005: ",sqrt(mean((finals[1].-finals[2]).^2)))
println("Final-field RMS dt .005/.0025: ",sqrt(mean((finals[2].-finals[3]).^2)))
c256, phi256 = initial_fields(256)
final,hist = run_simulation(c256,phi256;steps=1000,directory=joinpath(@__DIR__,"results-256"))
@test maximum(abs.(hist[:,2].-mean(c256))) < 1e-12
@test maximum(diff(hist[:,3])) < 1e-8
println("256² t=10: mass drift=$(maximum(abs.(hist[:,2].-mean(c256)))) energy=$(hist[1,3]) → $(hist[end,3]) range=$(extrema(final))")
