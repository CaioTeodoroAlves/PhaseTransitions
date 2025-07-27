module Dynamics

using ..Lattices
using Random

export DynamicsRule, RandomFlip, MajorityDynamics, HardcoreModel, IndependentResample, IsingDynamics
export step!, run_dynamics!
export random_config

# Abstract type for dynamics rules
abstract type DynamicsRule end

# Random flip dynamics - independently flip sites with given probability
struct RandomFlip <: DynamicsRule
    flip_probability::Float64
end

# Majority dynamics - each site updates based on majority of neighbors
struct MajorityDynamics <: DynamicsRule end

# Hardcore model - sites can only be active if no neighbors are active
struct HardcoreModel <: DynamicsRule end

# Ising dynamics (Glauber-style random updates)
struct IsingDynamics <: DynamicsRule
    beta::Float64  # Inverse temperature
    J::Int         # Coupling: +1 (ferro), -1 (antiferro)
end

# Independent resample rule (like site percolation)
struct IndependentResample <: DynamicsRule
    p::Float64
end

# Generic configuration type that works with any D-dimensional lattice
struct LatticeConfig{D}
    values::Dict{Any, Union{Bool, Int8}}
    lattice::RegularLattice{D}
end

# Constructor for LatticeConfig
function LatticeConfig(lattice::RegularLattice{D}, initial_values::Dict{Any, Union{Bool, Int8}}) where D
    LatticeConfig{D}(initial_values, lattice)
end

# Alternative constructor with different argument order
function LatticeConfig(initial_values::Dict{Any, Union{Bool, Int8}}, lattice::RegularLattice{D}) where D
    LatticeConfig{D}(initial_values, lattice)
end

# Helper function to get value at a site
function get_value(config::LatticeConfig{D}, site) where D
    return get(config.values, site, false)  # Default to false for missing sites
end

# Helper function to set value at a site
function set_value!(config::LatticeConfig{D}, site, value::Union{Bool, Int8}) where D
    config.values[site] = value
end

# Helper function to get all sites in the lattice
function all_sites(config::LatticeConfig{D}) where D
    return sites(config.lattice)
end

# Helper function to get random site
function random_site(config::LatticeConfig{D}) where D
    sites_list = collect(all_sites(config))
    return sites_list[rand(1:length(sites_list))]
end

"""
    step!(config, lattice, rule, update_fraction=1.0)

Apply one step of the given dynamics rule to the configuration, updating a fraction of sites at random.
"""
function step!(config::LatticeConfig{D}, rule::RandomFlip; update_fraction=1.0) where D
    sites_list = collect(all_sites(config))
    num_updates = round(Int, update_fraction * length(sites_list))
    for _ in 1:num_updates
        site = random_site(config)
        if rand() < rule.flip_probability
            current_value = get_value(config, site)
            set_value!(config, site, !current_value)
        end
    end
end

function step!(config::LatticeConfig{D}, rule::MajorityDynamics; update_fraction=1.0) where D
    sites_list = collect(all_sites(config))
    num_updates = round(Int, update_fraction * length(sites_list))
    for _ in 1:num_updates
        site = random_site(config)
        neighs = neighbors(config.lattice, site)
        if !isempty(neighs)
            occupied_neighbors = sum(get_value(config, neighbor) for neighbor in neighs)
            total_neighbors = length(neighs)
            # Update based on majority
            if occupied_neighbors > total_neighbors / 2
                set_value!(config, site, true)
            elseif occupied_neighbors < total_neighbors / 2
                set_value!(config, site, false)
            end
            # If exactly half, keep current state
        end
    end
end

function step!(config::LatticeConfig{D}, rule::HardcoreModel; update_fraction=1.0) where D
    sites_list = collect(all_sites(config))
    num_updates = round(Int, update_fraction * length(sites_list))
    for _ in 1:num_updates
        site = random_site(config)
        neighs = neighbors(config.lattice, site)
        if !isempty(neighs)
            has_occupied_neighbor = any(get_value(config, neighbor) for neighbor in neighs)
            if get_value(config, site) && has_occupied_neighbor
                set_value!(config, site, false)
            end
        end
    end
end

function step!(config::LatticeConfig{D}, rule::IsingDynamics; update_fraction=1.0) where D
    sites_list = collect(all_sites(config))
    num_updates = round(Int, update_fraction * length(sites_list))
    for _ in 1:num_updates
        site = random_site(config)
        s = get_value(config, site)
        neighs = neighbors(config.lattice, site)
        sum_neigh = sum(get_value(config, neighbor) for neighbor in neighs)
        dE = 2 * rule.J * s * sum_neigh
        if dE <= 0 || rand() < exp(-rule.beta * dE)
            set_value!(config, site, Int8(-s))
        end
    end
end

function step!(config::LatticeConfig{D}, rule::IndependentResample; update_fraction=1.0) where D
    sites_list = collect(all_sites(config))
    num_updates = round(Int, update_fraction * length(sites_list))
    for _ in 1:num_updates
        site = random_site(config)
        set_value!(config, site, rand() < rule.p)
    end
end

"""
    run_dynamics!(config, rule, steps; update_fraction=1.0)

Run the dynamics for a specified number of steps, updating a fraction of sites at random per step.
"""
function run_dynamics!(config::LatticeConfig, rule, steps; update_fraction=1.0)
    for _ in 1:steps
        step!(config, rule; update_fraction=update_fraction)
    end
end

# Helper to generate a random configuration for a given rule and lattice
function random_config(rule::MajorityDynamics, lattice::RegularLattice{D}) where D
    values = Dict{Any, Union{Bool, Int8}}()
    for site in sites(lattice)
        values[site] = rand(Bool)
    end
    return LatticeConfig(values, lattice)
end

function random_config(rule::RandomFlip, lattice::RegularLattice{D}) where D
    values = Dict{Any, Union{Bool, Int8}}()
    for site in sites(lattice)
        values[site] = rand(Bool)
    end
    return LatticeConfig(values, lattice)
end

function random_config(rule::HardcoreModel, lattice::RegularLattice{D}) where D
    values = Dict{Any, Union{Bool, Int8}}()
    for site in sites(lattice)
        values[site] = rand(Bool)
    end
    return LatticeConfig(values, lattice)
end

function random_config(rule::IndependentResample, lattice::RegularLattice{D}) where D
    values = Dict{Any, Union{Bool, Int8}}()
    for site in sites(lattice)
        values[site] = rand() < rule.p
    end
    return LatticeConfig(values, lattice)
end

function random_config(rule::IsingDynamics, lattice::RegularLattice{D}) where D
    values = Dict{Any, Union{Bool, Int8}}()
    for site in sites(lattice)
        values[site] = Int8(rand([-1, 1]))
    end
    return LatticeConfig(values, lattice)
end

# Backward compatibility: functions that accept lattice as separate argument
function step!(config::LatticeConfig{D}, lattice::RegularLattice{D}, rule; update_fraction=1.0) where D
    step!(config, rule; update_fraction=update_fraction)
end

function run_dynamics!(config::LatticeConfig{D}, lattice::RegularLattice{D}, rule, steps; update_fraction=1.0) where D
    run_dynamics!(config, rule, steps; update_fraction=update_fraction)
end

end # module 