/*
 * RunConfiguration.java
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

package dr.app.beast;

import dr.app.util.Arguments;
import dr.inference.mcmcmc.MCMCMCOptions;

/**
 * How the chains of a run are arranged: the number of independent chains, the parallel
 * tempering set run for each of them (a single chain when there are no heated chains) and
 * how many chains may compute at once.
 *
 * @author Andrew Rambaut
 */
public class RunConfiguration {

    public static final double DEFAULT_DELTA = 1.0;
    public static final int DEFAULT_SWAP_CHAIN_EVERY = 100;
    public static final int DEFAULT_CONVERGENCE_INTERVAL = 60;

    public RunConfiguration(int independentChainCount, double[] chainTemperatures, int swapChainsEvery,
                            MCMCMCOptions.SwapScheme swapScheme, int chainThreads) {
        if (independentChainCount < 1) {
            throw new IllegalArgumentException("At least one chain is needed");
        }
        if (chainTemperatures.length < 1 || chainTemperatures[0] != 1.0) {
            throw new IllegalArgumentException("The first chain must be cold (temperature 1.0)");
        }
        this.independentChainCount = independentChainCount;
        this.chainTemperatures = chainTemperatures;
        this.swapChainsEvery = swapChainsEvery;
        this.swapScheme = swapScheme;
        this.chainThreads = chainThreads;
    }

    /**
     * @param seconds the interval between convergence reports (0 for none)
     */
    public void setConvergenceInterval(int seconds) {
        this.convergenceInterval = seconds;
    }

    /**
     * @return the interval in seconds between convergence reports when there are several
     * independent chains (0 for none)
     */
    public int getConvergenceInterval() {
        return convergenceInterval;
    }

    /**
     * A single chain with no tempering.
     */
    public static RunConfiguration singleChain() {
        return new RunConfiguration(1, new double[]{1.0}, DEFAULT_SWAP_CHAIN_EVERY,
                MCMCMCOptions.SwapScheme.ORIGINAL_FLAVOR, -1);
    }

    /**
     * Reads the -chains, -chain_threads and -mc3_* options.
     *
     * @throws IllegalArgumentException if the options are inconsistent
     */
    public static RunConfiguration fromArguments(Arguments arguments) {

        int temperedChainCount = 1;
        if (arguments.hasOption("mc3_chains")) {
            temperedChainCount = arguments.getIntegerOption("mc3_chains");
        } else if (arguments.hasOption("mc3_temperatures")) {
            temperedChainCount = 1 + arguments.getRealArrayOption("mc3_temperatures").length;
        }

        double delta = DEFAULT_DELTA;
        if (arguments.hasOption("mc3_delta")) {
            if (arguments.hasOption("mc3_temperatures")) {
                throw new IllegalArgumentException("Either the -mc3_delta or the -mc3_temperatures option should be used, not both");
            }
            delta = arguments.getRealOption("mc3_delta");
        }

        double[] chainTemperatures = new double[temperedChainCount];
        chainTemperatures[0] = 1.0;
        if (arguments.hasOption("mc3_temperatures")) {
            double[] hotChainTemperatures = arguments.getRealArrayOption("mc3_temperatures");
            if (hotChainTemperatures.length != temperedChainCount - 1) {
                throw new IllegalArgumentException("The number of -mc3_temperatures should be one less than -mc3_chains");
            }
            System.arraycopy(hotChainTemperatures, 0, chainTemperatures, 1, temperedChainCount - 1);
        } else {
            for (int i = 1; i < temperedChainCount; i++) {
                chainTemperatures[i] = 1.0 / (1.0 + (delta * i));
            }
        }

        int swapChainsEvery = DEFAULT_SWAP_CHAIN_EVERY;
        if (arguments.hasOption("mc3_swap")) {
            swapChainsEvery = arguments.getIntegerOption("mc3_swap");
        }

        MCMCMCOptions.SwapScheme swapScheme = MCMCMCOptions.SwapScheme.ORIGINAL_FLAVOR;
        if (arguments.hasOption("mc3_scheme")) {
            swapScheme = MCMCMCOptions.SwapScheme.parse(arguments.getStringOption("mc3_scheme"));
        }

        int independentChainCount = 1;
        if (arguments.hasOption("chains")) {
            independentChainCount = arguments.getIntegerOption("chains");
        }

        int chainThreads = -1;
        if (arguments.hasOption("chain_threads")) {
            chainThreads = arguments.getIntegerOption("chain_threads");
        }

        RunConfiguration configuration = new RunConfiguration(independentChainCount, chainTemperatures,
                swapChainsEvery, swapScheme, chainThreads);
        if (arguments.hasOption("convergence_every")) {
            configuration.setConvergenceInterval(arguments.getIntegerOption("convergence_every"));
        }
        return configuration;
    }

    /**
     * @return true if this is a plain single-chain run
     */
    public boolean isSingleChain() {
        return getTotalChainCount() == 1;
    }

    /**
     * @return the number of chains in each parallel tempering set (1 means no tempering)
     */
    public int getTemperedChainCount() {
        return chainTemperatures.length;
    }

    /**
     * @return the number of independent chains (each with its tempering set if any)
     */
    public int getIndependentChainCount() {
        return independentChainCount;
    }

    /**
     * @return the number of chains in the whole run
     */
    public int getTotalChainCount() {
        return independentChainCount * getTemperedChainCount();
    }

    /**
     * @return the number of chains that may compute at once: as given, or by default as many
     * as there are processors
     */
    public int getChainThreads() {
        int total = getTotalChainCount();
        if (chainThreads > 0) {
            return Math.min(chainThreads, total);
        }
        return Math.min(total, Runtime.getRuntime().availableProcessors());
    }

    /**
     * @return the number of threads each computing chain should use for its likelihoods so
     * the run as a whole uses about one thread per processor
     */
    public int getThreadsPerChain() {
        return Math.max(1, Runtime.getRuntime().availableProcessors() / getChainThreads());
    }

    public MCMCMCOptions getMCMCMCOptions() {
        return new MCMCMCOptions(chainTemperatures.clone(), swapChainsEvery, swapScheme);
    }

    public double[] getChainTemperatures() {
        return chainTemperatures.clone();
    }

    private final int independentChainCount;
    private final double[] chainTemperatures;
    private final int swapChainsEvery;
    private final MCMCMCOptions.SwapScheme swapScheme;
    private final int chainThreads;
    private int convergenceInterval = DEFAULT_CONVERGENCE_INTERVAL;
}
