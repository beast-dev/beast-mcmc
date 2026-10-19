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
 * How the chains of a run are arranged: the parallel tempering set run for each cold chain
 * (a single chain when there are no heated chains) and how many chain threads may compute at
 * once.
 *
 * @author Andrew Rambaut
 */
public class RunConfiguration {

    public static final double DEFAULT_DELTA = 1.0;
    public static final int DEFAULT_SWAP_CHAIN_EVERY = 100;

    public RunConfiguration(double[] chainTemperatures, int swapChainsEvery,
                            MCMCMCOptions.SwapScheme swapScheme, int chainThreads) {
        if (chainTemperatures.length < 1 || chainTemperatures[0] != 1.0) {
            throw new IllegalArgumentException("The first chain must be cold (temperature 1.0)");
        }
        this.chainTemperatures = chainTemperatures;
        this.swapChainsEvery = swapChainsEvery;
        this.swapScheme = swapScheme;
        this.chainThreads = chainThreads;
    }

    /**
     * A single chain with no tempering.
     */
    public static RunConfiguration singleChain() {
        return new RunConfiguration(new double[]{1.0}, DEFAULT_SWAP_CHAIN_EVERY,
                MCMCMCOptions.SwapScheme.ORIGINAL_FLAVOR, -1);
    }

    /**
     * Reads the -mc3_* options.
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

        return new RunConfiguration(chainTemperatures, swapChainsEvery, swapScheme, -1);
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
     * @return the number of chains in the whole run
     */
    public int getTotalChainCount() {
        return getTemperedChainCount();
    }

    /**
     * @return the number of chains that may compute at once
     */
    public int getChainThreads() {
        int total = getTotalChainCount();
        if (chainThreads > 0) {
            return Math.min(chainThreads, total);
        }
        return total;
    }

    public MCMCMCOptions getMCMCMCOptions() {
        return new MCMCMCOptions(chainTemperatures.clone(), swapChainsEvery, swapScheme);
    }

    public double[] getChainTemperatures() {
        return chainTemperatures.clone();
    }

    private final double[] chainTemperatures;
    private final int swapChainsEvery;
    private final MCMCMCOptions.SwapScheme swapScheme;
    private final int chainThreads;
}
