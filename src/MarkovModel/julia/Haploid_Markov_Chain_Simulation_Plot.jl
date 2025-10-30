using Distributions
using LinearAlgebra
using Statistics 
using Plots # NEW: Package for plotting simulation results

# --- 1. Utility Functions ---

"""
    hamming_distance(g1::Int, g2::Int)

Calculates the Hamming distance between two binary genotypes (represented as integers).
This is the number of positions at which the two bits differ.
"""
function hamming_distance(g1::Int, g2::Int)::Int
    # XORing the two numbers highlights the differing bits (set to 1)
    # count_ones counts the number of set bits (1s)
    return count_ones(xor(g1, g2))
end

"""
    multinomial_pdf(N::Int, counts::Vector{Int}, probs::Vector{Float64})::Float64

Calculates the probability density function (PDF) of drawing a specific 'counts' 
vector (the next state) from the Multinomial distribution. (Used for calculating the T matrix).
"""
function multinomial_pdf(N::Int, counts::Vector{Int}, probs::Vector{Float64})::Float64
    return pdf(Multinomial(N, probs), counts)
end


# --- 2. Core Model Functions (Simulation Helper) ---

"""
    mutation_matrix(L::Int, mu::Float64)::Matrix{Float64}

Calculates the G x G matrix (M) where G = 2^L.
M[i, j] is the probability that genotype i mutates to genotype j.
This matrix defines the deterministic effect of mutation.
"""
function mutation_matrix(L::Int, mu::Float64)::Matrix{Float64}
    G = 2^L # Total number of genotypes
    M = zeros(G, G)

    for i in 0:(G - 1)
        for j in 0:(G - 1)
            d = hamming_distance(i, j)

            # P(i -> j) = mu^d * (1 - mu)^(L - d)
            M[i + 1, j + 1] = mu^d * (1.0 - mu)^(L - d)
        end
    end
    return M
end

"""
    simulate_generations(N::Int, L::Int, mu::Float64, T_gen::Int, X_initial::Vector{Int})

Runs the stochastic haploid population model for T_gen generations.

The transition from X_t (counts) to X_t+1 (counts) involves two steps:
1. Mutation: Calculate expected frequencies P_next = (X_t / N)' * M
2. Drift (Sampling): X_t+1 is sampled from Multinomial(N, P_next)
"""
function simulate_generations(N::Int, L::Int, mu::Float64, T_gen::Int, X_initial::Vector{Int})
    G = 2^L # Number of genotypes
    M = mutation_matrix(L, mu)
    
    # Store the count vectors and frequency vectors over time
    X_history = Vector{Vector{Int}}()
    P_history = Vector{Vector{Float64}}()

    X_current = X_initial
    
    # Check if initial state is valid
    if sum(X_initial) != N || length(X_initial) != G
        error("Initial state is invalid. Must sum to N=$N and have $G genotypes.")
    end

    push!(X_history, X_current)
    push!(P_history, X_current ./ N)

    for t in 1:T_gen
        # Step 1: Determine P_next_gen (Frequencies after Mutation and Random Mating)
        p_current = X_current ./ N
        
        # P_next_gen is the row vector of expected frequencies
        P_next_gen = p_current' * M 
        
        # We must convert the 1xG Adjoint to a standard Vector for Multinomial input
        P_next_vec = vec(P_next_gen) 

        # Step 2: Stochastic Sampling (Genetic Drift) via Multinomial
        # The next state X_next is a random sample of N individuals using P_next_vec probabilities.
        dist = Multinomial(N, P_next_vec)
        X_next = rand(dist)

        # Record and update
        push!(X_history, X_next)
        push!(P_history, X_next ./ N)
        X_current = X_next
    end

    return X_history, P_history
end


# --- 3. Example Execution (Simulation) ---

# Define Parameters: These are kept small for quick demonstration
const N = 50     # Population size (small N means high drift)
const L = 4       # Sequence length (4 genotypes: 00, 01, 10, 11)
const MU = 0.05    # Mutation rate per locus
const T_GEN = 100  # Number of generations to simulate

G = 2^L

# Initial State: 100% Genotype 0 (00)
# Genotype order: 00, 01, 10, 11 (index 1, 2, 3, 4)
X_initial = zeros(Int, G)
X_initial[1] = N # Start with N individuals of genotype 00

# Run the simulation
X_history, P_history = simulate_generations(N, L, MU, T_GEN, X_initial)

println("--- Stochastic Hybrid Model Simulation (N=$N, L=$L, mu=$MU) ---")

println("\nGenotype Mapping:")
for i in 0:(G - 1)
    # Using i+1 for 1-based Julia index
    println("Genotype Index $(i+1): $(string(i, base=2, pad=L))")
end

println("\nSimulation Results (First 5 and Last 5 Generations):")

# Print initial states
println("Gen 0 (Counts): $(X_history[1]) | Frequencies: $(round.(P_history[1], digits=3))")

# Print early generations
for t in 1:min(4, T_GEN)
    println("Gen $t (Counts): $(X_history[t+1]) | Frequencies: $(round.(P_history[t+1], digits=3))")
end

if T_GEN > 10
    println("...")
end

# Print late generations
for t in (T_GEN - min(5, T_GEN) + 1):T_GEN
    # t + 1 is the generation index (since history includes t=0)
    println("Gen $t (Counts): $(X_history[t+1]) | Frequencies: $(round.(P_history[t+1], digits=3))")
end

# --- 4. Plotting the Results ---
println("\n--- Generating Plot of Genotype Frequencies Over Time ---")

# 1. Prepare data for plotting
# P_history is Vector{Vector{Float64}}. We need a (T_GEN + 1) x G matrix 
# where each column is a genotype's frequency over time.
num_gens = T_GEN + 1
P_matrix = hcat(P_history...)' # Stack the vectors column-wise, then transpose

# Create labels for the legend
genotype_labels = [string(i, base=2, pad=L) for i in 0:(G - 1)]

# Plot the frequencies (columns of P_matrix) against the generation index (rows)
gr() # Use GR backend
p = plot(0:T_GEN, P_matrix, 
         label=reshape(genotype_labels, 1, :), # Reshape for legend labels
         xlabel="Generation", 
         ylabel="Genotype Frequency",
         title="Haploid Population Dynamics (N=$N, L=$L, μ=$MU)",
         legend=:right,
         linewidth=2,
         framestyle=:box,
         size=(800, 600)) # Set size for better display

# Display the plot
display(p)

println("\nFinal State (Gen $T_GEN):")
println("Counts: $(X_history[end])")
println("Frequencies: $(round.(P_history[end], digits=5))")

