using Random
using Statistics
using LinearAlgebra
using Printf
using Plots

"""
Moran model with ASSORTATIVE mating based on genetic similarity Q
Q = 1 - h/n where h is Hamming distance

Dispersal is implemented using a Kernel, directional dispersal is possible
"""
mutable struct MoranPopulationAssortative
    L::Int                          # Grid size (L x L)
    N::Int                          # Total number of individuals
    n::Int                          # Sequence length (number of bits)
    mu::Float64                     # Mutation rate per bit
    sigma::Float64                  # Dispersal distance
    Q_min::Float64                  # Minimum genetic similarity for mating (0 to 1)
    gamma::Float64                  # Directional bias strength
    theta0::Float64                 # Preferred direction in radians
    positions::Vector{Tuple{Int,Int}}
    sequences::Matrix{Int}
    time::Float64
    
    function MoranPopulationAssortative(;L=20, n=5, mu=0.01, sigma=3.0, Q_min=0.6, 
                                       gamma=0.0, theta0=0.0)
        N = L * L
        positions = [(i % L + 1, div(i, L) + 1) for i in 0:N-1]
        sequences = zeros(Int, N, n)
        time = 0.0
        new(L, N, n, mu, sigma, Q_min, gamma, theta0, positions, sequences, time)
    end
end

"""
Calculate Hamming distance between sequences of individuals i and j
"""
function hamming_distance(pop::MoranPopulationAssortative, i::Int, j::Int)
    return sum(pop.sequences[i, :] .!= pop.sequences[j, :])
end

"""
Calculate genetic similarity Q = 1 - h/n
Q = 1: identical (h=0)
Q = 0: maximally different (h=n)
"""
function genetic_similarity(pop::MoranPopulationAssortative, i::Int, j::Int)
    h = hamming_distance(pop, i, j)
    return 1.0 - h / pop.n
end

"""
Check if two individuals can mate
ASSORTATIVE: Only if Q ≥ Q_min (similar enough)
"""
function can_mate(pop::MoranPopulationAssortative, i::Int, j::Int; 
                  ignore_spatial=false)
    Q = genetic_similarity(pop, i, j)
    
    # Must be genetically similar enough
    if Q < pop.Q_min
        return false
    end
    
    # Optionally check spatial constraint
    if !ignore_spatial
        xi, yi = pop.positions[i]
        xj, yj = pop.positions[j]
        dist = sqrt((xi - xj)^2 + (yi - yj)^2)
        
        if dist > 3 * pop.sigma
            return false
        end
    end
    
    return true
end

"""
Mating probability between individuals i and j
Combines spatial distance and genetic similarity
"""
function mating_probability(pop::MoranPopulationAssortative, i::Int, j::Int)
    # Spatial component
    xi, yi = pop.positions[i]
    xj, yj = pop.positions[j]
    dx = xj - xi
    dy = yj - yi
    distance = sqrt(dx^2 + dy^2)
    
    if distance == 0
        angle = 0.0
    else
        angle = atan(dy, dx)
    end
    
    # Directional factor
    directional = 1.0 + pop.gamma * cos(angle - pop.theta0)
    eff_distance = distance / directional
    spatial_prob = exp(-eff_distance^2 / (2 * pop.sigma^2))
    
    # Genetic similarity constraint
    Q = genetic_similarity(pop, i, j)
    
    # ASSORTATIVE: Can only mate if similar enough (Q ≥ Q_min)
    if Q < pop.Q_min
        return 0.0  # Too different, cannot mate
    else
        return spatial_prob  # Similar enough, can mate
    end
end

"""
Build compatibility graph
Returns adjacency matrix where A[i,j] = 1 if i and j can mate
"""
function build_compatibility_graph(pop::MoranPopulationAssortative; 
                                   ignore_spatial=false)
    A = zeros(Int, pop.N, pop.N)
    
    for i in 1:pop.N
        for j in 1:pop.N
            if i != j && can_mate(pop, i, j, ignore_spatial=ignore_spatial)
                A[i, j] = 1
            end
        end
    end
    
    return A
end

"""
Find reproductive clusters using Union-Find algorithm
This identifies groups that can interbreed (species)
"""
function find_reproductive_clusters(pop::MoranPopulationAssortative; 
                                    ignore_spatial=false)
    # Build compatibility graph
    A = build_compatibility_graph(pop, ignore_spatial=ignore_spatial)
    
    # Union-Find data structure
    parent = collect(1:pop.N)
    rank = zeros(Int, pop.N)
    
    function find_root(x)
        if parent[x] != x
            parent[x] = find_root(parent[x])  # Path compression
        end
        return parent[x]
    end
    
    function union!(x, y)
        root_x = find_root(x)
        root_y = find_root(y)
        
        if root_x != root_y
            # Union by rank
            if rank[root_x] < rank[root_y]
                parent[root_x] = root_y
            elseif rank[root_x] > rank[root_y]
                parent[root_y] = root_x
            else
                parent[root_y] = root_x
                rank[root_x] += 1
            end
        end
    end
    
    # Build connected components
    for i in 1:pop.N
        for j in (i+1):pop.N
            if A[i, j] == 1
                union!(i, j)
            end
        end
    end
    
    # Identify clusters
    clusters = Dict{Int, Vector{Int}}()
    for i in 1:pop.N
        root = find_root(i)
        if !haskey(clusters, root)
            clusters[root] = Int[]
        end
        push!(clusters[root], i)
    end
    
    return collect(values(clusters))
end

"""
Count number of reproductive clusters (species)
"""
function count_species(pop::MoranPopulationAssortative; ignore_spatial=false)
    clusters = find_reproductive_clusters(pop, ignore_spatial=ignore_spatial)
    return length(clusters)
end

"""
Get detailed cluster information
"""
function analyze_clusters(pop::MoranPopulationAssortative; ignore_spatial=false)
    clusters = find_reproductive_clusters(pop, ignore_spatial=ignore_spatial)
    
    cluster_info = []
    
    for (idx, cluster) in enumerate(clusters)
        # Get phenotypes in this cluster
        phenotypes = [sum(pop.sequences[i, :]) for i in cluster]
        
        # Get sequences
        sequences_in_cluster = [pop.sequences[i, :] for i in cluster]
        
        # Calculate within-cluster genetic similarity
        Q_values = Float64[]
        for i in 1:(length(cluster)-1)
            for j in (i+1):length(cluster)
                Q = genetic_similarity(pop, cluster[i], cluster[j])
                push!(Q_values, Q)
            end
        end
        
        mean_Q = length(Q_values) > 0 ? mean(Q_values) : 1.0
        
        # Calculate cluster statistics
        info = Dict(
            "cluster_id" => idx,
            "size" => length(cluster),
            "individuals" => cluster,
            "mean_phenotype" => mean(phenotypes),
            "phenotype_range" => (minimum(phenotypes), maximum(phenotypes)),
            "genetic_diversity" => length(unique([join(s) for s in sequences_in_cluster])),
            "mean_similarity" => mean_Q
        )
        
        push!(cluster_info, info)
    end
    
    # Sort by size (largest first)
    sort!(cluster_info, by = x -> x["size"], rev=true)
    
    return cluster_info
end

"""
Calculate species diversity index (Simpson's)
"""
function species_diversity_index(cluster_sizes::Vector{Int})
    N_total = sum(cluster_sizes)
    if N_total == 0
        return 0.0
    end
    
    D = 1.0 - sum((n / N_total)^2 for n in cluster_sizes)
    return D
end

"""
Select mate for individual i
"""
function select_mate(pop::MoranPopulationAssortative, i::Int)
    probs = zeros(Float64, pop.N)
    
    for j in 1:pop.N
        probs[j] = mating_probability(pop, i, j)
    end
    
    prob_sum = sum(probs)
    
    
    # Normalize to probability distribution
    probs ./= prob_sum
    
    # Sample mate
    cumsum_probs = cumsum(probs)
    r = rand()
    
    for j in 1:pop.N
        if r <= cumsum_probs[j]
            return j
        end
    end
    
    return pop.N
end

"""
One Moran step: single birth-death event
"""
function moran_step!(pop::MoranPopulationAssortative)
    # Birth: select parent
    parent = rand(1:pop.N)
    mate = select_mate(pop, parent)
    
    # Create offspring via recombination
    offspring_sequence = zeros(Int, pop.n)
    for k in 1:pop.n
        if rand() < 0.5
            offspring_sequence[k] = pop.sequences[parent, k]
        else
            offspring_sequence[k] = pop.sequences[mate, k]
        end
    end
    
    # Apply mutations
    for k in 1:pop.n
        if rand() < pop.mu
            offspring_sequence[k] = 1 - offspring_sequence[k]
        end
    end
    
    # Death: select individual to replace
    dead = rand(1:pop.N)
    pop.sequences[dead, :] = offspring_sequence
    
    pop.time += 1.0
end

"""
One generation (N Moran steps)
"""
function moran_generation!(pop::MoranPopulationAssortative)
    for _ in 1:pop.N
        moran_step!(pop)
    end
end

"""
Calculate phenotype (sum of bits) for each individual
"""
function get_phenotypes(pop::MoranPopulationAssortative)
    return vec(sum(pop.sequences, dims=2))
end

"""
Get 2D spatial pattern of phenotypes
"""
function get_spatial_pattern(pop::MoranPopulationAssortative)
    phenotypes = get_phenotypes(pop)
    pattern = reshape(phenotypes, pop.L, pop.L)
    return pattern
end

"""
Get distribution of genetic similarities in population
"""
function similarity_distribution(pop::MoranPopulationAssortative)
    Q_values = Float64[]
    
    for i in 1:pop.N
        for j in (i+1):pop.N
            Q = genetic_similarity(pop, i, j)
            push!(Q_values, Q)
        end
    end
    
    return Q_values
end

"""
Print detailed cluster analysis
"""
function print_cluster_analysis(pop::MoranPopulationAssortative, 
                               generation::Int; ignore_spatial=false)
    println("\n" * "="^70)
    println("REPRODUCTIVE CLUSTER ANALYSIS - Generation $generation")
    println("ASSORTATIVE MATING (Q_min = $(pop.Q_min))")
    println("Time: $(pop.time / pop.N) generations")
    if ignore_spatial
        println("(Ignoring spatial constraints)")
    end
    println("="^70)
    
    cluster_info = analyze_clusters(pop, ignore_spatial=ignore_spatial)
    num_clusters = length(cluster_info)
    
    println("\nNumber of reproductive clusters (species): $num_clusters")
    
    # Diversity index
    cluster_sizes = [info["size"] for info in cluster_info]
    diversity_idx = species_diversity_index(cluster_sizes)
    @printf("Simpson's diversity index: %.3f\n", diversity_idx)
    
    # Overall genetic similarity distribution
    Q_values = similarity_distribution(pop)
    @printf("Population-wide genetic similarity: mean=%.3f, std=%.3f\n",
            mean(Q_values), std(Q_values))
    @printf("  Q ∈ [%.3f, %.3f]\n", minimum(Q_values), maximum(Q_values))
    
    println("\nCluster details:")
    println("-"^70)
    println("ID | Size | %Pop | Mean_P | P_Range | Diversity | Mean_Q")
    println("-"^70)
    
    for info in cluster_info
        @printf("%2d | %4d | %4.1f%% | %6.2f | [%d,%d] | %9d | %.3f\n",
                info["cluster_id"],
                info["size"],
                info["size"] / pop.N * 100,
                info["mean_phenotype"],
                info["phenotype_range"][1],
                info["phenotype_range"][2],
                info["genetic_diversity"],
                info["mean_similarity"])
    end
    
    # Spatial distribution visualization
    if num_clusters > 1 && num_clusters <= 10
        println("\nSpatial distribution of largest clusters:")
        pattern = fill(0, pop.L, pop.L)
        
        for (cluster_idx, info) in enumerate(cluster_info[1:min(5, num_clusters)])
            for ind in info["individuals"]
                x, y = pop.positions[ind]
                pattern[x, y] = cluster_idx
            end
        end
        
        display_size = min(15, pop.L)
        for i in 1:display_size
            for j in 1:display_size
                if pattern[i, j] == 0
                    print(". ")
                else
                    print("$(pattern[i, j]) ")
                end
            end
            println()
        end
        println("Legend: numbers = cluster ID, . = other clusters")
    end
    
    return cluster_info
end

"""
Track species dynamics over time
"""
function track_species_dynamics(;L=20, n=5, mu=0.01, sigma=3.0, Q_min=0.6,
                               gamma=0.0, theta0=0.0, generations=500,
                               sample_interval=50, ignore_spatial=false)
    
    println("="^70)
    println("TRACKING SPECIES DYNAMICS OVER TIME")
    println("ASSORTATIVE MATING: Q = 1 - h/n")
    println("="^70)
    println("\nParameters:")
    println("  Grid: $(L)x$(L), N=$(L*L)")
    println("  Sequence length: $n")
    println("  Mutation rate: $mu")
    println("  Dispersal σ: $sigma")
    println("  Min genetic similarity Q_min: $Q_min")
    println("    (corresponds to h_max = $(Int(round((1-Q_min)*n))))")
    println("  Generations: $generations")
    println("  Sampling interval: $sample_interval")
    
    pop = MoranPopulationAssortative(L=L, n=n, mu=mu, sigma=sigma, 
                                    Q_min=Q_min, gamma=gamma, theta0=theta0)
    
    # History tracking
    time_points = Int[]
    num_species_history = Int[]
    diversity_index_history = Float64[]
    mean_Q_history = Float64[]
    
    print_cluster_analysis(pop, 0, ignore_spatial=ignore_spatial)
    
    println("\n" * "="^70)
    println("Running simulation...")
    println("="^70)
    
    for gen in 1:generations
        moran_generation!(pop)
        
        if gen % sample_interval == 0
            cluster_info = analyze_clusters(pop, ignore_spatial=ignore_spatial)
            num_species = length(cluster_info)
            
            cluster_sizes = [info["size"] for info in cluster_info]
            diversity_idx = species_diversity_index(cluster_sizes)
            
            Q_values = similarity_distribution(pop)
            mean_Q = mean(Q_values)
            
            push!(time_points, gen)
            push!(num_species_history, num_species)
            push!(diversity_index_history, diversity_idx)
            push!(mean_Q_history, mean_Q)
            
            @printf("Gen %4d: %2d species, diversity=%.3f, mean_Q=%.3f\n",
                    gen, num_species, diversity_idx, mean_Q)
        end
    end
    
    print_cluster_analysis(pop, generations, ignore_spatial=ignore_spatial)
    
    # Summary
    println("\n" * "="^70)
    println("SPECIES DYNAMICS SUMMARY")
    println("="^70)
    
    println("\nEvolution over time:")
    println("Gen | Species | Diversity | Mean_Q")
    println("-"^45)
    for i in 1:length(time_points)
        @printf("%4d | %7d | %9.3f | %.3f\n",
                time_points[i],
                num_species_history[i],
                diversity_index_history[i],
                mean_Q_history[i])
    end
    
    println("\nFinal statistics:")
    @printf("  Initial species: %d\n", num_species_history[1])
    @printf("  Final species: %d\n", num_species_history[end])
    @printf("  Maximum species: %d\n", maximum(num_species_history))
    @printf("  Final diversity: %.3f\n", diversity_index_history[end])
    @printf("  Final mean Q: %.3f\n", mean_Q_history[end])
    
    return pop, time_points, num_species_history, diversity_index_history, mean_Q_history
end

"""
Compare speciation with different Q_min values
"""
function compare_Q_min_values(;L=15, n=5, mu=0.01, sigma=3.0, 
                              generations=400, sample_interval=100)
    
    println("\n" * "#"^70)
    println("# COMPARING SPECIATION ACROSS DIFFERENT Q_min VALUES")
    println("#"^70)
    
    # Test different Q_min values
    Q_min_values = [0.7, 0.8, 0.9, 0.95, 0.97]
    results = Dict()
    
    for Q_min in Q_min_values
        h_max_equiv = Int(round((1-Q_min)*n))
        println("\n" * ">"^70)
        println("> Q_min = $Q_min (equivalent to h_max = $h_max_equiv)")
        println(">"^70)
        
        pop, time_pts, num_species, diversity, mean_Q = track_species_dynamics(
            L=L, n=n, mu=mu, sigma=sigma, Q_min=Q_min,
            generations=generations, sample_interval=sample_interval,
            ignore_spatial=false
        )
        
        results[Q_min] = Dict(
            "final_species" => num_species[end],
            "max_species" => maximum(num_species),
            "final_diversity" => diversity[end],
            "final_mean_Q" => mean_Q[end]
        )
    end
    
    # Summary comparison
    println("\n" * "="^70)
    println("SPECIATION COMPARISON SUMMARY")
    println("="^70)
    println("\nQ_min | h_max | Final_Species | Max_Species | Diversity | Mean_Q")
    println("-"^70)
    
    for Q_min in Q_min_values
        h_max_equiv = Int(round((1-Q_min)*n))
        r = results[Q_min]
        @printf("%5.2f | %5d | %13d | %11d | %9.3f | %.3f\n",
                Q_min, h_max_equiv, r["final_species"], 
                r["max_species"], r["final_diversity"], r["final_mean_Q"])
    end
    
    return results
end


# ============================================================================
# SIMULATIONS
# ============================================================================



Random.seed!(123)


# Example 1: Basic simulation with Q_min = 0.96

pop1, time1, species1, diversity1, meanQ1 = track_species_dynamics(
    L=20, n=5, mu=0.01, sigma=3.0, Q_min=0.96,
    generations=400, sample_interval=100
)

# Example 2: Compare different Q_min values

results = compare_Q_min_values(
    L=15, n=5, mu=0.01, sigma=3.0,
    generations=400, sample_interval=100
)


