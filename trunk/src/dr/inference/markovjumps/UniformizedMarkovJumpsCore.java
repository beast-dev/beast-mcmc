package dr.inference.markovjumps;

/**
 * @author Marc Suchard
 * @author Vladimir Minin
 *
 * A class extension for implementing Markov chain-induced counting processes (markovjumps)
 * via uniformization in BEAST using BEAGLE
 * This work is supported by NSF grant 0856099
 *
 * Minin VN and Suchard MA (2008) Counting labeled transitions in continous-time Markov models of evolution.
 * Journal of Mathematical Biology, 56, 391-412.
 *
 * Rodrigue N, Philippe H and Lartillot N (2006) Uniformization for sampling realizations of Markov processes:
 * applications to Bayesian implementations of codon substitution models. Bioinformatics, 24, 56-62.
 *
 * Hobolth A and Stone E (2009) Simulation from endpoint-conditioned, continuous-time Markov chains on a finite
 * state space, with applications to molecular evolution. Annals of Applied Statistics, 3, 1204-1231.
 *
 */

public class UniformizedMarkovJumpsCore extends MarkovJumpsCore {

    public UniformizedMarkovJumpsCore(int stateCount) {
        super(stateCount);
    }

    // TODO Compute statistics of non-real-diagonalizable CTMC
    

}

/*

2.3. Uniformization
The final strategy that we consider permits sampling from X(t) through construction of an auxilliary stochastic process Y(t). Let ? = maxc Qc and define the process Y(t) by letting the state changes be determined by a discrete-time Markov process with transition matrix

(2.7)
Note that, by construction, we allow virtual state changes in which a jump occurs but the state does not change. Indeed, virtual state changes for state a are possible if Raa > 0. Next, let the epochs of state changes be determined by an independent Poisson process with rate ?. The stochastic process Y(t) is called a Markov chain subordinated to a Poisson process and is equivalent to the original continuous-time Markov chain X(t) as the following calculation shows:

(2.8)
This approach is commonly referred to as uniformization, and we adopt that language here. In what follows, we describe how uniformization can be used to construct an algorithm for exact sampling from X(t), conditional on the beginning and ending states.
It follows directly from (2.8) that the transition function of the Markov chain subordinated to a Poisson process is given by

Thus, the number of state changes N (including the virtual) for the conditional process that starts in X(0) = a and ends in X(T) = b is given by

(2.9)
Given the number of state changes N = n, the times t1, Щ, tn at which those state changes occur are uniformly distributed in the time interval [0, T]. Furthermore, the state changes X(t1), Щ, X(tn?1) are determined by a Markov chain with transition matrix R conditional on the beginning state X(0) = a and ending state X(tn) = b.
Putting these things together, we have the following algorithm for simulating a continuous-time Markov chain {X(t): 0 В t В T} conditional on the starting state X(0) = a and ending state X(T) = b.
Algorithm 5 (Uniformization)
Simulate the number of state changes n from the distribution (2.9).
If the number of state changes is 0, we are done: X(t) = a, 0 В t В T.
If the number of state changes is 1 and a = b, we are done: X(t) = a, 0 В t В T.
If the number of state changes is 1 and a ­ b simulate t1 uniformly random in [0, T], we are done: X(t) = a, t < t1, and X(t) = b, t1 В t В T.
When the number of state changes n is at least 2, simulate n independent uniform random numbers in [0, T] and sort the numbers in increasing order to obtain the times of state changes 0 < t1 < ссс< tn < T. Simulate X(t1), Щ, X(tn?1) from a discrete-time Markov chain with transition matrix R and conditional on starting state X(0) = a and ending state X(tn) = b. Determine which state changes are virtual and return the remaining changes and corresponding times of change.
Remark 6
In Step 1 above, we find the number of state changes n by simulating u from a Uniform(0, 1) distribution and letting n be the first time the cumulative sum of (2.9) exceeds u. When calculating the cumulative sum we need to raise R to powers 1 through n. These powers of R are stored because they are required in Step 5 of the algorithm. We use the eigenvalue decomposition (2.2) of Q to calculate Pab(t).
Remark 7
In Step 5 above we simulate X(ti), i = 1, Щ, n ? 1, from the discrete distribution with probability masses

 */
