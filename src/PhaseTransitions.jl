module PhaseTransitions

include("Lattices.jl")
using .Lattices
export AbstractLattice, RegularLattice, LatticeEdge, sites, neighbors, edges, point_to_int, int_to_point
export SquareLattice, TriangularLattice, CubicLattice, NearestNeighborLattice
export SquareNearestNeighborLattice, CubicNearestNeighborLattice

include("Percolation.jl")
using .Percolation
export SitePercResult, BondPercResult, run_site_percolation, run_bond_percolation, percolation_from_site_config, percolation_from_bond_config

include("Visualization.jl")
using .Visualization
export plot_clusters

include("Dynamics.jl")
using .Dynamics
export DynamicsRule, RandomFlip, MajorityDynamics, HardcoreModel, IndependentResample, IsingDynamics, step!, run_dynamics!, random_config
export LatticeConfig, get_value, set_value!, all_sites, random_site

end # module 