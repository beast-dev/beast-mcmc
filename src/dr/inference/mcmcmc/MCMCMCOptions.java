/*
 * MCMCMCOptions.java
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
import dr.inference.operators.OperatorSchedule;

/**
 * A class that brings together the auxillary information associated
 * with an ParallelMCMC analysis.
 *
 *
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @author Marc A. Suchard
 */
public class MCMCMCOptions {

    public MCMCMCOptions(final double[] temperatures, final int swapChainsEvery,
                         final SwapScheme swapScheme) {
        this.temperatures = temperatures;
        this.swapChainsEvery = swapChainsEvery;
        this.swapScheme = swapScheme;
    }

    public double[] getChainTemperatures() {
        return temperatures;
    }

    public int getSwapChainsEvery() {
        return swapChainsEvery;
    }

    public SwapScheme getSwapScheme() {  return swapScheme; }

    public enum SwapScheme {

        ORIGINAL_FLAVOR("original") {
            @Override
            public ParallelTempering factory(MarkovChain[] chains, OperatorSchedule[] schedules,
                                             MCMCMCOptions mcmcmcOptions) {
                return new ParallelTempering.OriginalFlavor(chains, schedules, mcmcmcOptions);
            }
        },
        STOCHASTIC_SINGLE("stochastic_single"){
            @Override
            public ParallelTempering factory(MarkovChain[] chains, OperatorSchedule[] schedules,
                                             MCMCMCOptions mcmcmcOptions) {
                return new ParallelTempering.StochasticSingleSwap(chains, schedules, mcmcmcOptions);
            }
        },
        STOCHASTIC_MULTIPLE("stochastic_multiple") {
            @Override
            public ParallelTempering factory(MarkovChain[] chains, OperatorSchedule[] schedules,
                                             MCMCMCOptions mcmcmcOptions) {
                return new ParallelTempering.StochasticMultipleSwap(chains, schedules, mcmcmcOptions);
            }
        },
        DETERMINISTIC_SINGLE("deterministic_single") {
            @Override
            public ParallelTempering factory(MarkovChain[] chains, OperatorSchedule[] schedules,
                                             MCMCMCOptions mcmcmcOptions) {
                return new ParallelTempering.DeterministicSingleSwap(chains, schedules, mcmcmcOptions);
            }
        },
        DETERMINISTIC_MULTIPLE("deterministic_multiple") {
            @Override
            public ParallelTempering factory(MarkovChain[] chains, OperatorSchedule[] schedules,
                                             MCMCMCOptions mcmcmcOptions) {
                return new ParallelTempering.DeterministicMultipleSwap(chains, schedules, mcmcmcOptions);
            }
        };

        SwapScheme(String name) {
            this.name = name;
        }

        public static SwapScheme parse(String schemeName) {
            for (SwapScheme scheme : SwapScheme.values()) {
                if (scheme.name.equalsIgnoreCase(schemeName)) {
                    return scheme;
                }
            }
            throw new RuntimeException("Unknown swap scheme '" + schemeName + "'");
        }

        private final String name;

        abstract public ParallelTempering factory(MarkovChain[] chains,
                                         OperatorSchedule[] schedules,
                                         MCMCMCOptions mcmcmcOptions);
    }

    private final double[] temperatures;
    private final int swapChainsEvery;
    private final SwapScheme swapScheme;
}
