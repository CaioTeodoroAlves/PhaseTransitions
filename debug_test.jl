using PhaseTransitions.Lattices
using PhaseTransitions.Percolation
using Random

println("=== Percolation Refactor Debug Test ===")

N = 10
p = 0.5
boundary = :periodic
seed = 42
Random.seed!(seed)

# Test 2D square nearest neighbor lattice
lat = SquareNearestNeighborLattice((N, N), boundary=boundary)

# Site percolation
site_result = run_site_percolation(lat, p)
println("Site percolation clusters: ", site_result.clusters)
println("Occupied sites (sum): ", count(site_result.occupied_sites))

# Bond percolation
bond_result = run_bond_percolation(lat, p)
println("Bond percolation clusters: ", bond_result.clusters)
println("Occupied bonds (sum): ", count(bond_result.occupied_edges))

# Test on 3D cubic lattice
lat3d = CubicNearestNeighborLattice((N, N, N), boundary=boundary)
site_result3d = run_site_percolation(lat3d, p)
println("3D Site percolation clusters: ", site_result3d.clusters)
println("3D Occupied sites (sum): ", count(site_result3d.occupied_sites))

println("\nDebug test completed!")