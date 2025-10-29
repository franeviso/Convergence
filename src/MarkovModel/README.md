Mutation random mating with infinite pop size

When you remove the population size $N$ from the Markov chain model, you are effectively taking the limit as **$N \to \infty$**. In this infinite-population limit, the random process (genetic drift) vanishes, and the model becomes entirely **deterministic**.

The system is no longer a stochastic Markov *Chain* but a deterministic set of **difference equations** governing the change in genotype frequencies. We can describe the entire process with a simple matrix multiplication.

I'll provide the Julia code for this deterministic model, where the state is now the vector of genotype frequencies.

### Deterministic Model ($\mathbf{N \to \infty}$)

1.  **State:** The state $P_t$ is a vector of genotype frequencies $P_t = (p_0, p_1, \dots, p_{2^L-1})$, where $p_i$ is the frequency of genotype $i$, and $\sum p_i = 1$.
2.  **Transition:** The transition is determined only by the mutation matrix $M$.
    $$\mathbf{P}_{t+1} = \mathbf{P}_t \mathbf{M}$$
    Where:
    * $\mathbf{P}_t$ is the row vector of frequencies at time $t$.
    * $\mathbf{M}$ is the $2^L \times 2^L$ mutation transition matrix, where $M_{i, j} = P(\text{genotype } i \to \text{genotype } j)$.

Here is the Julia implementation for this deterministic, frequency-based dynamics.


http://googleusercontent.com/immersive_entry_chip/0

### Summary of the Model

1.  **Matrix `M_transition`:** This acts as the "transition matrix." It captures the probability that any single individual's genetic sequence changes due to mutation from generation to generation.
2.  **Random Mating:** With $N=\infty$, random mating simply means the probability of sampling a parent of type $i$ is exactly $p_i$.
3.  **Deterministic Update:** The new frequency vector is found by multiplying the current frequency vector by the matrix $M$: $\mathbf{P}_{t+1} = \mathbf{P}_t \mathbf{M}$.

This deterministic model is often used as a good approximation when the population size ($N$) is very large and genetic drift is negligible compared to mutation.
