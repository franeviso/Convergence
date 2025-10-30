using Distributions
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

"""
    multinomial_coefficient(N::Int, counts::Vector{Int})

Calculates the multinomial coefficient N! / (n1! * n2! * ...).
Uses log space for numerical stability, though for small N standard factorials might suffice.
We use the Multinomial distribution PDF from the Distributions package, which is cleaner.
"""
function multinomial_pdf(N::Int, counts::Vector{Int}, probs::Vector{Float64})::Float64
    # The Multinomial distribution models the probability of observing 'counts'
    # given N trials and genotype probabilities 'probs'.
    return pdf(Multinomial(N, probs), counts)
end


# --- 2. Core Model Functions ---

"""
    mutation_matrix(L::Int, mu::Float64)::Matrix{Float64}

Calculates the G x G matrix (M) where G = 2^L.
M[i, j] is the probability that genotype i mutates to genotype j.
Assumes independent, symmetric mutation (0 <-> 1) at each of the L loci.
"""
function mutation_matrix(L::Int, mu::Float64)::Matrix{Float64}
    G = 2^L # Total number of genotypes
    M = zeros(G, G)

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
    generate_state_space(N::Int, G::Int)::Vector{Vector{Int}}

Generates all possible count vectors (states) (n1, n2, ..., nG) such that sum(ni) = N.
This function recursively computes the compositions of N into G parts.
"""
function generate_state_space(N::Int, G::Int)::Vector{Vector{Int}}
    states = Vector{Vector{Int}}()

    # Recursive helper function
    function generate_compositions(k::Int, current_sum::Int, current_comp::Vector{Int})
        if k == G # Base case: Last component
            last_comp = N - current_sum
            if last_comp >= 0
                push!(states, vcat(current_comp, last_comp))
            end
            return
        end

        # Recursive step
        for i in 0:(N - current_sum)
            generate_compositions(k + 1, current_sum + i, vcat(current_comp, i))
        end
    end

    if G == 0
        return Vector{Vector{Int}}()
    end
    generate_compositions(1, 0, Int[])
    return states
end

"""
    build_transition_matrix(N::Int, L::Int, mu::Float64)

Constructs the full transition matrix M for the Markov chain.
T[X_i, X_j] is the probability of transitioning from state X_i to state X_j.
Returns the row-stochastic transition matrix and the state space vector.
"""
# REMOVED the incorrect type hint '::Matrix{Float64}' which caused the error
function build_transition_matrix(N::Int, L::Int, mu::Float64)
    G = 2^L # Number of genotypes (e.g., L=2 -> G=4)

    # Step 1: Generate the Genotype Mutation Matrix (M)
    M = mutation_matrix(L, mu)

    # Step 2: Generate the State Space (S)
    S = generate_state_space(N, G)
    num_states = length(S)
    println("Population Size N=$N, Sequence Length L=$L. Total Genotypes G=$G.")
    println("Total Markov Chain States: $num_states")

    # Step 3: Initialize the Transition Matrix (T)
    T = zeros(Float64, num_states, num_states)

    # Step 4: Calculate Transition Probabilities
    for i in 1:num_states
        X_current = S[i] # Current state: vector of counts (n1, n2, ...)
        
        # Calculate the marginal probability P(g_j) for each genotype in the next generation.
        # This combines random mating (p_i) and mutation (M).
        
        # Current genotype frequencies p = n/N (Random Mating part)
        p_current = X_current ./ N

        # The probability that a random offspring is genotype j, P(g_j), is:
        # P(g_j) = sum over i of [ P(parent is g_i) * P(g_i mutates to g_j) ]
        # This is equivalent to the matrix multiplication: p_current * M
        P_next_gen = p_current' * M
        
        # Safety check: ensure probabilities sum to 1 (within tolerance)
        if !isapprox(sum(P_next_gen), 1.0)
            @warn "Probabilities do not sum to 1.0 for state $i: $(sum(P_next_gen))"
        end

        # Now calculate the transition probability to every possible future state X_j
        for j in 1:num_states
            X_next = S[j] # Next state: vector of counts (n'1, n'2, ...)
            
            # P(X_current -> X_next) is the Multinomial probability
            # of drawing X_next (counts) given N trials and probabilities P_next_gen
            T[j, i] = multinomial_pdf(N, X_next, vec(P_next_gen))
        end
    end
    
    # The transition matrix T is constructed such that T[to, from] (column-stochastic)
    # We transpose it to be T[from, to] (row-stochastic) for standard Markov chain notation.
    T_row_stochastic = T'

    # Final check: all rows should sum to 1
    for i in 1:num_states
        if !isapprox(sum(T_row_stochastic[i, :]), 1.0)
            @warn "Row $i does not sum to 1.0: $(sum(T_row_stochastic[i, :]))"
        end
    end

    return T_row_stochastic, S
end


# --- 3. Example Execution ---

# Define Parameters: These MUST be small for the code to run in a reasonable time.
# N=3, L=2 results in G=4 genotypes and 20 states. This is a good illustrative size.
const N = 10      # Population size
const L = 5     # Sequence length (00, 01, 10, 11)
const MU = 0.01  # Mutation rate per locus (mu)

# Build the Transition Matrix
T_matrix, State_Space = build_transition_matrix(N, L, MU)

println("\n--- Genotype Mapping (L=$L) ---")
println("Genotype Index (0-based) -> Binary Sequence")
for i in 0:(2^L - 1)
    println("$i -> $(string(i, base=2, pad=L))")
end

println("\n--- State Space Mapping (Index -> Counts) ---")
for (i, state) in enumerate(State_Space)
    println("State $(i) (n00, n01, n10, n11) = $state")
end

println("\n--- Transition Matrix (T) Dimensions: $(size(T_matrix)) ---")
println("T[i, j] is P(State i -> State j)")
# Displaying the matrix for small N and L
println("T_matrix (partial view):")
display(T_matrix[1:min(5, size(T_matrix, 1)), 1:min(5, size(T_matrix, 2))])
# Print the full matrix only if it's very small
if size(T_matrix, 1) <= 10
    println("\nFull Transition Matrix:")
    display(T_matrix)
end

