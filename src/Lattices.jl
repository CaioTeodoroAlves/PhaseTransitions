module Lattices

import Base: size

# Exports
export AbstractLattice, SquareLattice
export sites, neighbors, point_to_int, int_to_point

# Abstract type for all lattice types
abstract type AbstractLattice end

# Concrete implementation of a square lattice
struct SquareLattice <: AbstractLattice
    N::Int
    boundary::Symbol  # :free or :periodic
    function SquareLattice(N::Int; boundary::Symbol=:free)
        new(N, boundary)
    end
end

size(lat::SquareLattice) = (lat.N, lat.N)
sites(lat::SquareLattice) = ((i, j) for i in 1:lat.N, j in 1:lat.N)

const NEIGHBOR_DIRECTIONS = ((-1, 0), (1, 0), (0, -1), (0, 1))

function neighbors(lat::SquareLattice, p::Tuple{Int, Int})
    i, j = p
    N = lat.N
    neighs = NTuple{2, Int}[]
    for (di, dj) in NEIGHBOR_DIRECTIONS
        ni, nj = i + di, j + dj
        if lat.boundary == :periodic
            ni = Int((ni - 1) % N + 1)
            nj = Int((nj - 1) % N + 1)
            # Robust fix: ensure indices are always in 1:N
            if ni == 0
                ni = N
            end
            if nj == 0
                nj = N
            end
            if !(1 <= ni <= N && 1 <= nj <= N)
                @warn "neighbors: Out-of-bounds after wrapping: ($ni, $nj) from ($i, $j) with boundary periodic"
            end
            push!(neighs, (ni, nj))
        elseif lat.boundary == :free
            if 1 <= ni <= N && 1 <= nj <= N
                push!(neighs, (ni, nj))
            end
        else
            error("Unknown boundary condition: $(lat.boundary)")
        end
    end
    return neighs
end

function point_to_int(lat::SquareLattice, p::Tuple{Int, Int})
    i, j = p
    return (j - 1) * lat.N + i
end

function int_to_point(lat::SquareLattice, n::Int)
    i = (n - 1) % lat.N + 1
    j = (n - 1) ÷ lat.N + 1
    return (i, j)
end

end # module 