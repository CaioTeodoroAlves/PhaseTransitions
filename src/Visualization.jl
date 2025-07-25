module Visualization

using ..Lattices
using ..Percolation
using Plots
using Images
using DataStructures
using ColorSchemes
using Random

export plot_clusters, save_cluster_animation

const MAX_PALETTE_SIZE = 1000
const _fixed_palette_cache = Dict{Any, Vector{RGB{Float64}}}()

function plot_clusters(result::SitePercResult; 
                       color_scheme=:random, 
                       palette=nothing,
                       saturation=0.85, 
                       value=0.95,
                       unoccupied_color=RGB(1,1,1),
                       color_seed=nothing)
    """
    Plot clusters with various color schemes.
    - color_scheme: :random (default), :hsv, :size_ordered, or palette
    - color_seed: if provided, use this integer to seed the RNG for color assignment (for reproducible colors)
    """
    lat = result.lattice
    clusters = result.clusters
    occupied = result.occupied_sites
    root_colors = Dict{Int, RGB{Float64}}()
    
    # Prepare palette if provided
    grad = palette !== nothing ? cgrad(palette) : nothing
    n_colors = palette !== nothing ? length(grad.colors) : 0

    # Set a fixed random seed for consistent coloring if requested
    rng = color_seed === nothing ? Random.GLOBAL_RNG : MersenneTwister(color_seed)

    dims = size(lat)
    if length(dims) != 2
        error("plot_clusters currently only supports 2D lattices.")
    end

    if color_scheme == :size_ordered
        # 1. Collect all cluster roots and their member sites
        cluster_sites = Dict{Int, Vector{Tuple{Int,Int}}}()
        for (i, j) in Iterators.product(1:dims[1], 1:dims[2])
            if occupied[i, j]
                p_int = point_to_int(lat, (i, j))
                root = find_root!(clusters, p_int)
                if !haskey(cluster_sites, root)
                    cluster_sites[root] = Vector{Tuple{Int,Int}}()
                end
                push!(cluster_sites[root], (i, j))
            end
        end
        # 2. Sort clusters by size (desc), then by top-left-most site
        sorted_roots = sort(collect(keys(cluster_sites)); lt=(r1, r2) -> begin
            s1, s2 = cluster_sites[r1], cluster_sites[r2]
            if length(s1) != length(s2)
                return length(s2) < length(s1)  # descending size
            else
                # break ties by top-left-most site
                min1 = findmin(s1)[1]
                min2 = findmin(s2)[1]
                return min1 < min2
            end
        end)
        # 3. Assign colors in order from a fixed-length palette
        n_clusters = length(sorted_roots)
        if grad !== nothing
            n_palette = length(grad.colors)
            color_list = Vector{RGB{Float64}}(undef, n_clusters)
            for idx in 1:n_clusters
                pal_idx = 1 + mod(idx-1, n_palette)
                color_list[idx] = grad.colors[pal_idx]
            end
        else
            cache_key = color_seed
            if haskey(_fixed_palette_cache, cache_key)
                palette_vec = _fixed_palette_cache[cache_key]
            else
                palette_vec = [rand(rng, RGB{Float64}) for _ in 1:MAX_PALETTE_SIZE]
                _fixed_palette_cache[cache_key] = palette_vec
            end
            # Dynamically grow the palette if needed
            if n_clusters > length(palette_vec)
                append!(palette_vec, [rand(rng, RGB{Float64}) for _ in 1:(n_clusters - length(palette_vec))])
                _fixed_palette_cache[cache_key] = palette_vec
            end
            color_list = palette_vec[1:n_clusters]
        end
        for (idx, root) in enumerate(sorted_roots)
            root_colors[root] = color_list[mod1(idx, length(color_list))]
        end
    else
        # Assign colors to cluster roots (original logic)
        for (i, j) in Iterators.product(1:dims[1], 1:dims[2])
            if occupied[i, j]
                p_int = point_to_int(lat, (i, j))
                root = find_root!(clusters, p_int)
                if !haskey(root_colors, root)
                    h = hash(root)
                    if grad !== nothing
                        # Use a float in [0,1] for palette lookup
                        colorval = mod1(h, 10^9) / 10^9
                        root_colors[root] = grad[colorval]
                    elseif color_scheme == :hsv
                        hue = h % 360
                        root_colors[root] = HSV(hue, saturation, value)
                    else
                        root_colors[root] = rand(rng, RGB{Float64})
                    end
                end
            end
        end
    end

    img = fill(unoccupied_color, dims) # User-specified background
    for (i, j) in Iterators.product(1:dims[1], 1:dims[2])
        if occupied[i, j]
            p_int = point_to_int(lat, (i, j))
            root = find_root!(clusters, p_int)
            if haskey(root_colors, root)
                img[i, j] = root_colors[root]
            end
        end
    end
    
    p = plot(img, aspect_ratio=:equal, axis=false, ticks=false)
    plot!(p, title="Cluster Plot")
    return p
end

function save_cluster_animation!(config, lat::AbstractLattice, rule, steps::Int; filename="cluster_animation.gif", interval=1, palette=nothing, color_scheme=:random, saturation=0.85, value=0.95, unoccupied_color=RGB(1,1,1))
    # config: BitMatrix (Boolean) or Matrix{Int8} (Ising)
    # rule: any dynamics rule with step! defined
    # steps: total number of frames
    # interval: number of dynamics steps between frames
    # filename: output GIF file
    anim = Plots.Animation()
    for t in 1:steps
        # Extract clusters for current config
        if typeof(config) <: BitMatrix
            perc = percolation_from_config(lat, config)
        elseif typeof(config) <: Matrix{Int8}
            perc = percolation_from_config(lat, config)
        else
            error("Unsupported config type: $(typeof(config))")
        end
        p = plot_clusters(perc; palette=palette, color_scheme=color_scheme, saturation=saturation, value=value, unoccupied_color=unoccupied_color)
        Plots.frame(anim, p)
        # Advance the dynamics
        for _ in 1:interval
            step!(config, lat, rule)
        end
    end
    gif(anim, filename, fps=15)
    return nothing
end

end # module 