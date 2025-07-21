module PhaseTransitions

include("Lattices.jl")
using .Lattices
export AbstractLattice, SquareLattice, sites, neighbors, point_to_int, int_to_point

include("Percolation.jl")
using .Percolation
export PercResult, run_site_percolation, percolation_from_config

include("Visualization.jl")
using .Visualization
export plot_clusters

include("Dynamics.jl")
using .Dynamics
export DynamicsRule, RandomFlip, MajorityDynamics, HardcoreModel, IndependentResample, IsingDynamics, step!, run_dynamics!, random_config

end # module 