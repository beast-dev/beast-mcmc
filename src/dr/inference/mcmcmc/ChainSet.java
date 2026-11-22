/*
 * ChainSet.java
 *
 * Copyright © 2002-2024 the BEAST Development Team
 * http://beast.community/about
 *
 * This file is part of BEAST.
 * See the NOTICE file distributed with this work for additional
 * information regarding copyright ownership and licensing.
 *
 * BEAST is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 *  BEAST is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with BEAST; if not, write to the
 * Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
 * Boston, MA  02110-1301  USA
 *
 */

package dr.inference.mcmcmc;

import dr.inference.markovchain.MarkovChain;
import dr.inference.mcmc.MCMC;

import java.util.concurrent.Semaphore;

/**
 * A group of MCMC chains that sample one posterior and are driven together by a
 * MultiChainRunner: a single chain, or a cold chain with its parallel tempering partners.
 * Sets are independent of each other; chains within a set may synchronise.
 *
 * @author Andrew Rambaut
 */
public interface ChainSet {

    /**
     * @return the index of this set among the sets being run
     */
    int getIndex();

    /**
     * Prepares the chains to run. Each chain gets its own thread and runs only while holding
     * one of the given permits, which bounds how many chains compute at once.
     *
     * @param permits shared by all the sets in a run
     */
    void start(Semaphore permits);

    /**
     * Runs the next unit of work (for a tempered set, the segment between two swaps). Called
     * repeatedly by the runner on the set's driver thread.
     *
     * @return true if there is more to run
     */
    boolean step();

    /**
     * Closes the loggers and reports. Only called when the set ran to the end (or was stopped).
     */
    void finish();

    /**
     * Releases the chain threads. Always called, also after a failure.
     */
    void shutdown();

    /**
     * Requests that all the chains in the set stop. Can be called from any thread.
     */
    void pleaseStop();

    /**
     * @return the MCMC currently sampling the target (temperature 1) distribution
     */
    MCMC getColdMCMC();

    /**
     * @return the Markov chain currently sampling the target distribution
     */
    MarkovChain getColdChain();

    /**
     * @return the number of states run so far
     */
    long getCurrentState();

    /**
     * @return a one-line status for progress reports
     */
    String getStatus();
}
