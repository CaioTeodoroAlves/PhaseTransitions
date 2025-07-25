# To run this from the Julia REPL:
# 1. cd to the PhaseTransitions directory
# 2. Start Julia: `julia`
# 3. Enter pkg mode by typing `]`
# 4. Activate the environment: `activate .`
# 5. Test the package: `test` or `include("test/runtests.jl")`

using PhaseTransitions
using Test
using DataStructures

@testset "PhaseTransitions.jl" begin

    println("--- Testing PhaseTransitions.jl ---")

    # 1. Test Lattice creation
    N = 10
    lat = SquareLattice(N)
    @test size(lat) == (N, N)
    println("✓ Lattice creation successful.")

    # 2. Test coordinate conversion
    p = (3, 4)
    n = point_to_int(lat, p)
    p_back = int_to_point(lat, n)
    @test p == p_back
    println("✓ Coordinate conversion successful.")

    # 3. Test percolation run
    p_crit = 0.5927
    result = run_site_percolation(lat, p_crit)
    @test typeof(result) == SitePercResult
    @test length(result.clusters.parents) == N * N
    num_occupied = sum(result.occupied_sites)
    final_clusters_small = num_groups(result.clusters)
    println("✓ Percolation run successful.")
    println("  - Lattice size: $N x $N")
    println("  - Occupation probability p = $p_crit")
    println("  - Number of occupied sites: $num_occupied")
    println("  - Final number of clusters: $final_clusters_small")

    # 4. Test percolation with larger N and check cluster count
    println("\n--- Testing with larger lattice ---")
    N_large = 100
    lat_large = SquareLattice(N_large)
    p_large = 0.5
    result_large = run_site_percolation(lat_large, p_large)

    num_occupied_large = sum(result_large.occupied_sites)
    final_clusters_large = num_groups(result_large.clusters)

    println("✓ Large percolation run successful.")
    println("  - Lattice size: $N_large x $N_large")
    println("  - Occupation probability p = $p_large")
    println("  - Number of occupied sites: $num_occupied_large")
    println("  - Final number of clusters: $final_clusters_large")

    @test final_clusters_large <= N_large * N_large
    @test final_clusters_large >= (N_large * N_large - num_occupied_large)

    println("\n--- All tests passed! ---")
    
end

@testset "SquareLattice neighbors (free)" begin
    lat = SquareLattice(5, boundary=:free)
    # Corner
    @test sort(neighbors(lat, (1,1))) == [(1,2), (2,1)]
    # Edge
    @test sort(neighbors(lat, (1,3))) == [(1,2), (1,4), (2,3)]
    # Center
    @test sort(neighbors(lat, (3,3))) == [(2,3), (3,2), (3,4), (4,3)]
end

@testset "SquareLattice neighbors (periodic)" begin
    lat = SquareLattice(5, boundary=:periodic)
    # Corner
    @test sort(neighbors(lat, (1,1))) == [(1,2), (2,1), (1,5), (5,1)]
    # Edge
    @test sort(neighbors(lat, (1,3))) == [(1,2), (1,4), (2,3), (5,3)]
    # Center
    @test sort(neighbors(lat, (3,3))) == [(2,3), (3,2), (3,4), (4,3)]
end

@testset "Percolation and Dynamics with boundaries" begin
    for boundary in (:free, :periodic)
        lat = SquareLattice(10, boundary=boundary)
        # Percolation
        result = run_site_percolation(lat, 0.5)
        @test size(result.occupied_sites) == (10,10)
        # Dynamics
        config = BitMatrix(rand(Bool, 10, 10))
        step!(config, lat, MajorityDynamics())
        step!(config, lat, HardcoreModel())
        step!(config, lat, IndependentResample(0.5))
        # Ising
        spins = Int8.(rand([-1,1], 10, 10))
        step!(spins, lat, IsingDynamics(0.5, 1))
        # Cluster extraction
        result2 = percolation_from_config(lat, config)
        @test size(result2.occupied_sites) == (10,10)
    end
end 