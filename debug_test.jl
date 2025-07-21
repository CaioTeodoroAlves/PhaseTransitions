using PhaseTransitions
using Plots
using Random

println("=== Spin Dynamics Video Demo ===")

N = 500
steps = 100
seed = 42
boundary = :periodic

# Choose dynamics rule: MajorityDynamics or IsingDynamics(beta, J)
dynamics_rule = IsingDynamics(0.44, -1)

# Initialize random configuration appropriate for the rule
debug_seed = seed  # for reproducibility
Random.seed!(debug_seed)
config = random_config(dynamics_rule, N)
lat = SquareLattice(N, boundary=boundary)
update_fraction = 0.1

# Prepare animation
anim = @animate for t in 1:steps
    step!(config, lat, dynamics_rule)
    # For visualization, map spins to Bool (e.g., up spins as true)
    perc_result = percolation_from_config(lat, config isa BitMatrix ? config : config .== 1)
    p = plot_clusters(perc_result, color_scheme=:size_ordered, palette=:viridis)
    plot!(p, title="Step $t", legend=false)
end

gif(anim, "dynamics_demo.gif", fps=3)
println("\nDemo completed! Animation saved as 'dynamics_demo.gif'")