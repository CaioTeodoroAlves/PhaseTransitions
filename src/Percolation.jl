module Percolation

using ..Lattices
using DataStructures

export PercResult, run_site_percolation, percolation_from_config

struct PercResult
    lattice::SquareLattice
    occupied_sites::BitMatrix
    clusters::IntDisjointSets
end

function run_site_percolation(lat::SquareLattice, p::Float64)
    occupied = rand(size(lat)...) .< p
    clusters = IntDisjointSets(lat.N * lat.N)
    N = lat.N

    for j in 1:N, i in 1:N # Iterate through each site
        if !occupied[i, j]
            continue
        end
        site_int = point_to_int(lat, (i, j))
        for (ni, nj) in neighbors(lat, (i, j))
            if occupied[ni, nj]
                neighbor_int = point_to_int(lat, (ni, nj))
                union!(clusters, site_int, neighbor_int)
            end
        end
    end
    return PercResult(lat, occupied, clusters)
end

"""
    percolation_from_config(lat::SquareLattice, config::BitMatrix)

Create a PercResult from a Boolean configuration, treating true as occupied sites.
"""
function percolation_from_config(lat::SquareLattice, config::BitMatrix)
    clusters = IntDisjointSets(lat.N * lat.N)
    N = lat.N
    for j in 1:N, i in 1:N
        if !config[i, j]
            continue
        end
        site_int = point_to_int(lat, (i, j))
        for (ni, nj) in neighbors(lat, (i, j))
            if !(1 <= ni <= N && 1 <= nj <= N)
                @warn "Out-of-bounds neighbor index: ($ni, $nj) for site ($i, $j) with boundary $(lat.boundary)"
                continue
            end
            if config[ni, nj]
                neighbor_int = point_to_int(lat, (ni, nj))
                union!(clusters, site_int, neighbor_int)
            end
        end
    end
    return PercResult(lat, config, clusters)
end

"""
    percolation_from_config(lat::SquareLattice, spins::Matrix{Int8})

Create a PercResult from an Ising spin configuration, treating +1 as occupied sites.
"""
function percolation_from_config(lat::SquareLattice, spins::Matrix{Int8})
    occupied = BitMatrix(spins .== 1)
    return percolation_from_config(lat, occupied)
end

end # module 