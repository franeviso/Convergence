Model explicitly **coupling** the deterministic effect of **mutation** (via the matrix $M$) with the stochastic effect of **genetic drift** (via the Multinomial distribution) to calculate the transition probabilities between states.

### Coupling

1.  **State Definition:** The Markov chain state $X$ is defined by the **count vector** $\mathbf{X} = (n_0, n_1, \dots, n_G)$, where $n_i$ is the number of individuals with genotype $i$, and $\sum n_i = N$.
2.  **Expected Frequencies ($P_{next}$):** Given the current state $\mathbf{X}_{current}$, the current frequencies are $\mathbf{P}_{current} = \mathbf{X}_{current} / N$. The *expected* frequency of genotype $j$ in the next generation, factoring in mutation but *not* drift, is:
    $$\mathbf{P}_{next} = \mathbf{P}_{current} \mathbf{M}$$
3.  **Stochastic Transition (Multinomial):** The actual next state $\mathbf{X}_{next}$ is not deterministic. It's determined by $N$ independent draws from the population, where the probability of drawing genotype $j$ is $P_{next, j}$. This random sampling process is described by the **Multinomial distribution**:
    $$P(\mathbf{X}_{current} \to \mathbf{X}_{next}) = \text{Multinomial}(\mathbf{X}_{next} \mid N, \mathbf{P}_{next})$$




