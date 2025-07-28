module Lattices

import Base: size

# Exports
export AbstractLattice, RegularLattice, LatticeEdge
export sites, neighbors, edges, point_to_int, int_to_point
export SQUARE_NEIGHBORHOOD, TRIANGULAR_NEIGHBORHOOD, SquareLattice, TriangularLattice
export CubicLattice
export NearestNeighborLattice

# Abstract type for all lattice types
abstract type AbstractLattice end

# General parametric lattice type
struct RegularLattice{D} <: AbstractLattice
    size::NTuple{D,Int}
    neighborhood::Vector{NTuple{D,Int}}  # relative offsets
    boundary::Symbol  # :free or :periodic
    function RegularLattice(size::NTuple{D,Int}, neighborhood::Vector{NTuple{D,Int}}; boundary::Symbol=:free) where D
        new{D}(size, neighborhood, boundary)
    end
end

# Edge representation
struct LatticeEdge{D}
    site1::NTuple{D,Int}
    site2::NTuple{D,Int}
end

size(lat::RegularLattice{D}) where D = lat.size

function sites(lat::RegularLattice{D}) where D
    ranges = ntuple(i -> 1:lat.size[i], D)
    return Iterators.product(ranges...)
end

function neighbors(lat::RegularLattice{D}, p::NTuple{D,Int}) where D
    neighs = NTuple{D,Int}[]
    for offset in lat.neighborhood
        neighbor = ntuple(i -> p[i] + offset[i], D)
        if lat.boundary == :periodic
            wrapped = ntuple(i -> mod(neighbor[i] - 1, lat.size[i]) + 1, D)
            push!(neighs, wrapped)
        elseif lat.boundary == :free
            inbounds = true
            for i in 1:D
                if neighbor[i] < 1 || neighbor[i] > lat.size[i]
                    inbounds = false
                    @warn "Neighbor $neighbor has out-of-bounds coordinate at dimension $i: $(neighbor[i]) not in [1, $(lat.size[i])]"
                    break
                end
            end
            if inbounds
                push!(neighs, neighbor)
            end
        else
            error("Unknown boundary condition: $(lat.boundary)")
        end
    end
    return neighs
end

function edges(lat::RegularLattice{D}) where D
    edge_set = Set{Tuple{NTuple{D,Int}, NTuple{D,Int}}}()
    for site in sites(lat)
        s = Tuple(site)
        for neighbor in neighbors(lat, s)
            edge = s < neighbor ? (s, neighbor) : (neighbor, s)
            push!(edge_set, edge)
        end
    end
    return [LatticeEdge{D}(e[1], e[2]) for e in edge_set]
end

function point_to_int(lat::RegularLattice{D}, p::NTuple{D,Int}) where D
    # Check bounds
    for d in 1:D
        if p[d] < 1 || p[d] > lat.size[d]
            error("Point $p is out of bounds for lattice of size $(lat.size)")
        end
    end
    idx = 1
    for d in 1:D
        stride = d == 1 ? 1 : prod(lat.size[1:d-1])
        idx += (p[d] - 1) * stride
    end
    return idx
end

function int_to_point(lat::RegularLattice{D}, n::Int) where D
    idx = n - 1
    coords = Int[]
    for d in 1:D
        stride = d == 1 ? 1 : prod(lat.size[1:d-1]...)
        coord = (idx ÷ stride) + 1
        push!(coords, coord)
        idx = idx % stride
    end
    return tuple(coords...)
end

const SQUARE_NEIGHBORHOOD = [
    ( 1,  0),
    (-1,  0),
    ( 0,  1),
    ( 0, -1)
]

"""
    SquareLattice(size::Tuple{Int,Int}; boundary=:free)

Construct a 2D square lattice of given size and boundary condition.
"""
function SquareLattice(size::Tuple{Int,Int}; boundary::Symbol=:free)
    RegularLattice(size, SQUARE_NEIGHBORHOOD; boundary=boundary)
end

const TRIANGULAR_NEIGHBORHOOD = [
    ( 1,  0),
    (-1,  0),
    ( 0,  1),
    ( 0, -1),
    ( 1, -1),
    (-1,  1)
]

"""
    TriangularLattice(size::Tuple{Int,Int}; boundary=:free)

Construct a 2D triangular lattice of given size and boundary condition.
"""
function TriangularLattice(size::Tuple{Int,Int}; boundary::Symbol=:free)
    RegularLattice(size, TRIANGULAR_NEIGHBORHOOD; boundary=boundary)
end

const CUBIC_NEIGHBORHOOD = [
    ( 1,  0,  0),
    (-1,  0,  0),
    ( 0,  1,  0),
    ( 0, -1,  0),
    ( 0,  0,  1),
    ( 0,  0, -1)
]

"""
    CubicLattice(size::Tuple{Int,Int,Int}; boundary=:free)

Construct a 3D cubic lattice of given size and boundary condition.
"""
function CubicLattice(size::Tuple{Int,Int,Int}; boundary::Symbol=:free)
    RegularLattice(size, CUBIC_NEIGHBORHOOD; boundary=boundary)
end

"""
    NearestNeighborLattice(size::NTuple{D,Int}; boundary=:free)

Construct a d-dimensional nearest-neighbor lattice of given size and boundary condition.
"""
function NearestNeighborLattice(size::NTuple{D,Int}; boundary::Symbol=:free) where D
    neighborhood = NTuple{D,Int}[]
    for d in 1:D
        vplus = ntuple(i -> i == d ? 1 : 0, D)
        vminus = ntuple(i -> i == d ? -1 : 0, D)
        push!(neighborhood, vplus)
        push!(neighborhood, vminus)
    end
    RegularLattice(size, neighborhood; boundary=boundary)
end

end # module 