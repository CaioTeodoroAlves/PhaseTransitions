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

"""
    step!(config, lattice, rule, update_fraction=1.0)

Apply one step of the given dynamics rule to the configuration, updating a fraction of sites at random.
"""
function step!(config::BitMatrix, lattice::RegularLattice{2}, rule::RandomFlip; update_fraction=1.0)
    dims = size(lattice)
    num_updates = round(Int, update_fraction * prod(dims))
    for _ in 1:num_updates
        i = rand(1:dims[1])
        j = rand(1:dims[2])
        if rand() < rule.flip_probability
            config[i, j] = !config[i, j]
        end
    end
end

function step!(config::BitMatrix, lattice::RegularLattice{2}, rule::MajorityDynamics; update_fraction=1.0)
    dims = size(lattice)
    num_updates = round(Int, update_fraction * prod(dims))
    for _ in 1:num_updates
        i = rand(1:dims[1])
        j = rand(1:dims[2])
        neighs = neighbors(lattice, (i, j))
        if !isempty(neighs)
            occupied_neighbors = sum(config[ni, nj] for (ni, nj) in neighs)
            total_neighbors = length(neighs)
            # Update based on majority
            if occupied_neighbors > total_neighbors / 2
                config[i, j] = true
            elseif occupied_neighbors < total_neighbors / 2
                config[i, j] = false
            end
            # If exactly half, keep current state
        end
    end
end

function step!(config::BitMatrix, lattice::RegularLattice{2}, rule::HardcoreModel; update_fraction=1.0)
    dims = size(lattice)
    num_updates = round(Int, update_fraction * prod(dims))
    for _ in 1:num_updates
        i = rand(1:dims[1])
        j = rand(1:dims[2])
        neighs = neighbors(lattice, (i, j))
        if !isempty(neighs)
            has_occupied_neighbor = any(config[ni, nj] for (ni, nj) in neighs)
            if config[i, j] && has_occupied_neighbor
                config[i, j] = false
            end
        end
    end
end

function step!(spins::Matrix{Int8}, lattice::RegularLattice{2}, rule::IsingDynamics; update_fraction=1.0)
    dims = size(lattice)
    num_updates = round(Int, update_fraction * prod(dims))
    for _ in 1:num_updates
        i = rand(1:dims[1])
        j = rand(1:dims[2])
        s = spins[i, j]
        neighs = neighbors(lattice, (i, j))
        sum_neigh = sum(spins[ni, nj] for (ni, nj) in neighs)
        dE = 2 * rule.J * s * sum_neigh
        if dE <= 0 || rand() < exp(-rule.beta * dE)
            spins[i, j] = -s
        end
    end
end

function step!(config::BitMatrix, lattice::RegularLattice{2}, rule::IndependentResample; update_fraction=1.0)
    dims = size(lattice)
    num_updates = round(Int, update_fraction * prod(dims))
    for _ in 1:num_updates
        i = rand(1:dims[1])
        j = rand(1:dims[2])
        config[i, j] = rand() < rule.p
    end
end

"""
    run_dynamics!(config, lattice, rule, steps; update_fraction=1.0)

Run the dynamics for a specified number of steps, updating a fraction of sites at random per step.
"""
function run_dynamics!(config, lattice, rule, steps; update_fraction=1.0)
    for _ in 1:steps
        step!(config, lattice, rule; update_fraction=update_fraction)
    end
end

# Helper to generate a random configuration for a given rule and size
function random_config(rule::MajorityDynamics, lattice::RegularLattice{2})
    BitMatrix(rand(Bool, size(lattice)...))
end

function random_config(rule::RandomFlip, lattice::RegularLattice{2})
    BitMatrix(rand(Bool, size(lattice)...))
end

function random_config(rule::HardcoreModel, lattice::RegularLattice{2})
    BitMatrix(rand(Bool, size(lattice)...))
end

function random_config(rule::IndependentResample, lattice::RegularLattice{2})
    BitMatrix(rand(Bool, size(lattice)...))
end

function random_config(rule::IsingDynamics, lattice::RegularLattice{2})
    Int8.(rand([-1, 1], size(lattice)...))
end

end # module 