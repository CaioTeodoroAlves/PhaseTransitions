module Percolation

using ..Lattices
using DataStructures

export SitePercResult, BondPercResult, run_site_percolation, run_bond_percolation, percolation_from_site_config, percolation_from_bond_config

struct SitePercResult{L<:Lattices.AbstractLattice, D}
    lattice::L
    occupied_sites::BitArray{D}
    clusters::IntDisjointSets
end

struct BondPercResult{L<:Lattices.AbstractLattice, D}
    lattice::L
    occupied_edges::BitVector
    edge_list::Vector{Lattices.LatticeEdge{D}}
    clusters::IntDisjointSets
end

function run_site_percolation(lat::Lattices.RegularLattice{D}, p::Float64) where D
    occ = rand(size(lat)...) .< p
    occ_bits = BitArray(occ)
    n_sites = prod(size(lat))
    clusters = IntDisjointSets(n_sites)
    for site in sites(lat)
        if !occ_bits[site...]
            continue
        end
        site_int = point_to_int(lat, Tuple(site))
        for neighbor in neighbors(lat, Tuple(site))
            if occ_bits[neighbor...]
                neighbor_int = point_to_int(lat, neighbor)
                union!(clusters, site_int, neighbor_int)
            end
        end
    end
    return SitePercResult(lat, occ_bits, clusters)
end

function run_bond_percolation(lat::Lattices.RegularLattice{D}, p::Float64) where D
    edge_list = edges(lat)
    n_edges = length(edge_list)
    occ = rand(n_edges) .< p
    clusters = IntDisjointSets(prod(size(lat)))
    for (idx, edge) in enumerate(edge_list)
        if !occ[idx]
            continue
        end
        site1_int = point_to_int(lat, edge.site1)
        site2_int = point_to_int(lat, edge.site2)
        union!(clusters, site1_int, site2_int)
    end
    return BondPercResult{typeof(lat), D}(lat, BitVector(occ), edge_list, clusters)
end

function percolation_from_site_config(lat::Lattices.RegularLattice{D}, config::BitArray{D}) where D
    n_sites = prod(size(lat))
    clusters = IntDisjointSets(n_sites)
    for site in sites(lat)
        if !config[site...]
            continue
        end
        site_int = point_to_int(lat, Tuple(site))
        for neighbor in neighbors(lat, Tuple(site))
            if config[neighbor...]
                neighbor_int = point_to_int(lat, neighbor)
                union!(clusters, site_int, neighbor_int)
            end
        end
    end
    return SitePercResult(lat, config, clusters)
end

function percolation_from_site_config(lat::Lattices.RegularLattice{D}, config::Dict) where D
    # Convert Dict-based config to BitArray
    occ_bits = falses(size(lat)...)
    for site in sites(lat)
        if haskey(config, site) && config[site]
            occ_bits[site...] = true
        end
    end
    return percolation_from_site_config(lat, occ_bits)
end

function percolation_from_bond_config(lat::Lattices.RegularLattice{D}, bond_config::BitVector) where D
    edge_list = edges(lat)
    clusters = IntDisjointSets(prod(size(lat)))
    for (idx, edge) in enumerate(edge_list)
        if !bond_config[idx]
            continue
        end
        site1_int = point_to_int(lat, edge.site1)
        site2_int = point_to_int(lat, edge.site2)
        union!(clusters, site1_int, site2_int)
    end
    return BondPercResult{typeof(lat), D}(lat, bond_config, edge_list, clusters)
end

end # module 