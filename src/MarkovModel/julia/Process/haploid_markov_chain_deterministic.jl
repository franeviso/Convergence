using LinearAlgebra

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

# --- 2. Core Model Function ---

"""
    mutation_matrix(L::Int, mu::Float64)::Matrix{Float64}

Calculates the G x G matrix (M) where G = 2^L.
M[i, j] is the probability that genotype i mutates to genotype j.
Assumes independent, symmetric mutation (0 <-> 1) at each of the L loci.
M is the transition matrix for a single genome.
"""
function mutation_matrix(L::Int, mu::Float64)::Matrix{Float64}
    G = 2^L # Total number of genotypes
    M = zeros(G, G)

    # i and j represent the 0-indexed genotypes (0 to G-1)
    for i in 0:(G - 1)
        for j in 0:(G - 1)
            # Calculate Hamming distance between genotype i and j
            d = hamming_distance(i, j)

            # Probability that d loci mutate and L-d loci do not mutate
            # P(i -> j) = mu^d * (1 - mu)^(L - d)
            M[i + 1, j + 1] = mu^d * (1.0 - mu)^(L - d)
        end
    end
    return M
end

"""
    next_generation(P_current::Vector{Float64}, M::Matrix{Float64})::Vector{Float64}

Calculates the genotype frequencies in the next generation (P_t+1)
based on the current frequencies (P_t) and the mutation matrix (M).

Transition Rule (Deterministic): P_t+1' = P_t' * M
(Where P_t is a row vector, M is the genotype transition matrix)
"""
function next_generation(P_current::Vector{Float64}, M::Matrix{Float64})::Vector{Float64}
    # P_current' * M performs the required row vector times matrix multiplication.
    # The result is a 1xG Adjoint vector, so we convert it back to a standard Vector{Float64}.
    P_next = (P_current' * M)'
    
    # Safety check
    if !isapprox(sum(P_next), 1.0)
        error("Frequencies do not sum to 1.0: $(sum(P_next))")
    end

    return P_next
end


# --- 3. Example Simulation ---

# Define Parameters
const L = 3      # Sequence length (8 genotypes: 000 to 111)
const MU = 0.05  # Mutation rate per locus per generation
const T_GEN = 10 # Number of generations to simulate

# Step 1: Genotype Mapping
G = 2^L
println("--- Genotype Mapping (L=$L) ---")
for i in 0:(G - 1)
    println("Genotype $i: $(string(i, base=2, pad=L))")
end
println("---------------------------------")

# Step 2: Initialize Frequencies (e.g., population starts 100% genotype 00...0)
P_initial = zeros(G)
P_initial[1] = 1.0 # Genotype 0 (00...0) has frequency 1.0
P_current = P_initial

# Step 3: Calculate the Mutation Matrix
M_transition = mutation_matrix(L, MU)

println("Initial Frequencies (Gen 0): $(P_initial)")
println("Simulating $T_GEN generations...")

# Step 4: Run Simulation
for t in 1:T_GEN
    P_next = next_generation(P_current, M_transition)
    
    # Optional: Print every few generations
    if t % 2 == 0 || t == T_GEN
        println("Gen $t Frequencies: $(round.(P_next, digits=4))")
    end
    
    global P_current = P_next
end

# Find the expected equilibrium (Stationary Distribution)
# In this simple, non-selective, symmetric mutation model, the equilibrium 
# distribution is uniform (p_i = 1/G).
P_eq = ones(G) ./ G
println("\nExpected Equilibrium (Uniform): $(round.(P_eq, digits=4))")
println("Distance from Equilibrium after $T_GEN generations: $(norm(P_current - P_eq))")

