/*
 * ParallelTempering.java
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
import dr.inference.mcmc.MCMCCriterion;
import dr.inference.operators.AdaptableMCMCOperator;
import dr.inference.operators.MCMCOperator;
import dr.inference.operators.OperatorSchedule;
import dr.math.MathUtils;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

/**
 * @author Marc A. Suchard
 * @author Guy Baele
 *
 * literature to look at:
 *   Saifuddin Syed and Alex Bouchard-Cote (see extensions to `New` below)
 *   Adapative parallel tempering for BEAST 2 (https://github.com/nicfel/CoupledMCMC)
 */

public interface ParallelTempering {

    int swapChainTemperatures(int coldChain);

    String getReport();

    // TODO implement various swapping approaches

    class IndexPair {
        int index1, index2;

        IndexPair(int index1, int index2) {
            this.index1 = index1;
            this.index2 = index2;
        }
    }

    abstract class Base implements ParallelTempering {

        final MarkovChain[] chains;
        final OperatorSchedule[] schedules;
        final ParallelTemperingStatistics statistics;

        Base(MarkovChain[] chains, OperatorSchedule[] schedules, MCMCMCOptions options) {
            this.chains = chains;
            this.schedules = schedules;
            this.statistics = new ParallelTemperingStatistics(options);
        }

        abstract List<IndexPair> getPairsToSwap();

        public int swapChainTemperatures(int coldChain) {

            List<IndexPair> pairs = getPairsToSwap();

            for (IndexPair pair : pairs) {

                boolean swap = scoreSwap(pair);
                if (swap) {
                    swapTemperatures(pair);
                    swapOperatorSchedules(pair);
                    coldChain = getNewColdChain(pair, coldChain);
                }
            }

            return coldChain;
        }

        public String getReport() {
            return statistics.getReport();
        }

        private int getNewColdChain(IndexPair pair, int coldChain) {
            if (pair.index1 == coldChain) {
                coldChain = pair.index2;
            } else if (pair.index2 == coldChain) {
                coldChain = pair.index1;
            }
            return coldChain;
        }

        private boolean scoreSwap(IndexPair pair) {

            double score1 = chains[pair.index1].getCurrentScore();
            MCMCCriterion acceptor1 = ((MCMCCriterion) chains[pair.index1].getAcceptor());
            double temperature1 = acceptor1.getTemperature();
            int rank1 = acceptor1.getRank();

            double score2 = chains[pair.index2].getCurrentScore();
            MCMCCriterion acceptor2 = ((MCMCCriterion) chains[pair.index2].getAcceptor());
            double temperature2 = acceptor2.getTemperature();
            int rank2 = acceptor2.getRank();

            double logRatio = ((score2 - score1) * temperature1) + ((score1 - score2) * temperature2);
            boolean success = (Math.log(MathUtils.nextDouble()) < logRatio);

            statistics.recordStatistics(pair.index1, pair.index2,
                    rank1, rank2,
                    temperature1, temperature2,
                    logRatio, success);

            return success;
        }

        private void swapTemperatures(IndexPair pair) {

            MCMCCriterion acceptor1 = ((MCMCCriterion) chains[pair.index1].getAcceptor());
            double temperature1 = acceptor1.getTemperature();

            MCMCCriterion acceptor2 = ((MCMCCriterion) chains[pair.index2].getAcceptor());
            double temperature2 = acceptor2.getTemperature();

            acceptor1.setTemperature(temperature2);
            acceptor2.setTemperature(temperature1);

            int rank1 = acceptor1.getRank();
            int rank2 = acceptor2.getRank();
            acceptor1.setRank(rank2);
            acceptor2.setRank(rank1);
        }

        private void swapOperatorSchedules(IndexPair pair) {

            OperatorSchedule schedule1 = schedules[pair.index1];
            OperatorSchedule schedule2 = schedules[pair.index2];

            for (int i = 0; i < schedule1.getOperatorCount(); i++) {
                MCMCOperator operator1 = schedule1.getOperator(i);
                MCMCOperator operator2 = schedule2.getOperator(i);

                long tmp = operator1.getAcceptCount();
                operator1.setAcceptCount(operator2.getAcceptCount());
                operator2.setAcceptCount(tmp);

                tmp = operator1.getRejectCount();
                operator1.setRejectCount(operator2.getRejectCount());
                operator2.setRejectCount(tmp);

                double tmp2 = operator1.getSumDeviation();
                operator1.setSumDeviation(operator2.getSumDeviation());
                operator2.setSumDeviation(tmp2);

                if (operator1 instanceof AdaptableMCMCOperator) {
                    tmp2 = ((AdaptableMCMCOperator) operator1).getAdaptableParameter();
                    ((AdaptableMCMCOperator) operator1).setAdaptableParameter(((AdaptableMCMCOperator) operator2).getAdaptableParameter());
                    ((AdaptableMCMCOperator) operator2).setAdaptableParameter(tmp2);
                }
            }
        }
    }

    class OriginalFlavor extends Base {

        public OriginalFlavor(MarkovChain[] chains, OperatorSchedule[] schedules, MCMCMCOptions options) {
            super(chains, schedules, options);
        }

        @Override
        List<IndexPair> getPairsToSwap() {
            int index1 = MathUtils.nextInt(chains.length);
            int index2 = MathUtils.nextInt(chains.length);
            while (index1 == index2) {
                index2 = MathUtils.nextInt(chains.length);
            }

            return Collections.singletonList(new IndexPair(index1, index2));
        }
    }

    abstract class New extends Base {

        private final ParitySelector paritySelector;
        private final HalfIndexSelector indexSelector;
        private final int K;
        private final int halfK;
        private boolean evenStep;

        public New(MarkovChain[] chains, OperatorSchedule[] schedules, MCMCMCOptions options,
                   ParitySelector directionSelector,
                   HalfIndexSelector indexSelector) {
            super(chains, schedules, options);
            this.paritySelector = directionSelector;
            this.indexSelector = indexSelector;

            int length = options.getChainTemperatures().length;
            this.K = length;
            this.halfK = length / 2;
            this.mapRankToChain = new int[K];

            this.evenStep = MathUtils.nextBoolean();
        }

        private static boolean isOdd(int i) {
            return (i & 1) == 1;
        }

        private final int[] mapRankToChain;

        @Override
        List<IndexPair> getPairsToSwap() {

            for (int i = 0; i < K; ++i) {
                MCMCCriterion acceptor = ((MCMCCriterion) chains[i].getAcceptor());
                int rank = acceptor.getRank();
                mapRankToChain[rank] = i;
            }

            evenStep = paritySelector.next(evenStep);
            int[] halfIndices = indexSelector.indices(halfK);

            List<IndexPair> pairs = new ArrayList<>();
            for (int halfIndex : halfIndices) {
                int rank1 = 2 * halfIndex + (evenStep ? 0 : 1);
                int rank2 = rank1 + 1;
                if (rank2 >= K) {
                    rank2 = 0;
                }

                pairs.add(new IndexPair(mapRankToChain[rank1], mapRankToChain[rank2]));
            }

            return pairs;
        }
    }

    class DeterministicSingleSwap extends New {
        public DeterministicSingleSwap(MarkovChain[] chains, OperatorSchedule[] schedules, MCMCMCOptions options) {
            super(chains, schedules, options, ParitySelector.DETERMINISTIC, HalfIndexSelector.ONE);
        }
    }

    class DeterministicMultipleSwap extends New {
        public DeterministicMultipleSwap(MarkovChain[] chains, OperatorSchedule[] schedules, MCMCMCOptions options) {
            super(chains, schedules, options, ParitySelector.DETERMINISTIC, HalfIndexSelector.ALL);
        }
    }

    class StochasticSingleSwap extends New {
        public StochasticSingleSwap(MarkovChain[] chains, OperatorSchedule[] schedules, MCMCMCOptions options) {
            super(chains, schedules, options, ParitySelector.STOCHASTIC, HalfIndexSelector.ONE);
        }
    }

    class StochasticMultipleSwap extends New {
        public StochasticMultipleSwap(MarkovChain[] chains, OperatorSchedule[] schedules, MCMCMCOptions options) {
            super(chains, schedules, options, ParitySelector.STOCHASTIC, HalfIndexSelector.ALL);
        }
    }

    enum ParitySelector {
        STOCHASTIC {
            @Override
            boolean next(boolean last) {
                return MathUtils.nextBoolean();
            }
        },
        DETERMINISTIC {
            @Override
            boolean next(boolean last) {
                return !last;
            }
        };

        abstract boolean next(boolean last);
    }

    enum HalfIndexSelector {
        ONE {
            @Override
            int[] indices(int size) {
                return new int[]{ MathUtils.nextInt(size) };
            }
        },
        ALL {
            @Override
            int[] indices(int size) {
                if (all == null) {
                    all = new int[size];
                    for (int i = 0; i < size; ++i) {
                        all[i] = i;
                    }
                }
                return all;
            }
        };

        abstract int[] indices(int size);

        private static int[] all = null;
    }
}
