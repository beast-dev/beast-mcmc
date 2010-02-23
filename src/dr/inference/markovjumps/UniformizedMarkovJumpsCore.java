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
