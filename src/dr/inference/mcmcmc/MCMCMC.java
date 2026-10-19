/*
 * MCMCMC.java
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

import dr.inference.loggers.LogFormatter;
import dr.inference.loggers.Logger;
import dr.inference.loggers.MCLogger;
import dr.inference.markovchain.MarkovChain;
import dr.inference.markovchain.MarkovChainListener;
import dr.inference.mcmc.MCMC;
import dr.inference.mcmc.MCMCCriterion;
import dr.inference.mcmc.MCMCOptions;
import dr.inference.model.Likelihood;
import dr.inference.model.Model;
import dr.inference.operators.OperatorAnalysisPrinter;
import dr.inference.operators.OperatorSchedule;
import dr.util.NumberFormatter;

import java.util.Collections;
import java.util.List;

/**
 * An MCMC analysis that estimates parameters of a probabilistic model.
 *
 * @author Andrew Rambaut
 */
public class MCMCMC implements Runnable {

    public final static boolean DEBUG = false;

    public MCMCMC(MCMC[] mcmcs, MCMCMCOptions mcmcmcOptions) {

        this.mcmcmcOptions = mcmcmcOptions;

        if (mcmcmcOptions.getChainTemperatures()[0] != 1.0) {
            throw new RuntimeException("The first chain in the array should be cold (temperature = 1.0)");
        }

        coldChain = 0;

        this.mcmcOptions = mcmcs[coldChain].getOptions();

        // Get all the loggers out of all the chains. We will only use the
        // loggers of the cold chain, but we need to swap the formatters around
        // so that which every chain is cold always writes to the same destination.
        mcLoggers = new MCLogger[mcmcs.length][];
        for (int i = 0; i < mcmcs.length; i++) {
            Logger[] loggers = mcmcs[i].getLoggers();
            mcLoggers[i] = new MCLogger[loggers.length];
            for (int j = 0; j < loggers.length; j++) {
                mcLoggers[i][j] = (MCLogger) loggers[j];
            }
            if (mcLoggers[i] == null) {
                throw new RuntimeException("There are no loggers in the MCMC chains.");
            }
        }

        // Get all the operator schedules. The tuning values of these must be swapped
        // around as the temperatures are swapped.
        schedules = new OperatorSchedule[mcmcs.length];
        for (int i = 0; i < schedules.length; i++) {
            schedules[i] = mcmcs[i].getOperatorSchedule();
        }

        chains = new MarkovChain[mcmcs.length];

        chains[0] = mcmcs[0].getMarkovChain();
        for (int i = 1; i < chains.length; i++) {
            chains[i] = mcmcs[i].getMarkovChain();
            MCMCCriterion acceptor = ((MCMCCriterion) chains[i].getAcceptor());
            acceptor.setTemperature(mcmcmcOptions.getChainTemperatures()[i]);
            acceptor.setRank(i);
        }

        scheme = mcmcmcOptions.getSwapScheme().factory(chains, schedules, mcmcmcOptions);
    }

    private static final boolean DEBUG_IN_SERIES = false;

    public void run() {
        currentState = 0;

        timer.start();

//        if (isPreBurninNeeded()) {
//            long preBurnin = mcmcOptions.getAdaptationDelay();
//            if (preBurnin > 0) {
//                MarkovChainListener burninListener = new BurninListener(preBurnin);
//
//                chains[coldChain].addMarkovChainListener(burninListener);
//                runChains(preBurnin, true);
//                chains[coldChain].removeMarkovChainListener(burninListener);
//                resetChains();
//            }
//        }

        MCLogger[] coldChainLoggers = mcLoggers[coldChain];
        List<LogFormatter>[] logFormatters = new List[coldChainLoggers.length];

        for (int i = 0; i < coldChainLoggers.length; i++) {
            // Start the logging for the cold chain
            coldChainLoggers[i].startLogging();

            // Now get the formatters (destinations) for the cold chains coldChainLoggers
            logFormatters[i] = coldChainLoggers[i].getFormatters();
        }

        // Set the other chains to have null log formatters...
        for (int j = 0; j < mcLoggers.length; j++) {
            if (j != coldChain) {
                for (int i = 0; i < mcLoggers[j].length; i++) {
                    mcLoggers[j][i].setFormatters(Collections.EMPTY_LIST);
                }
            }
        }

        chains[coldChain].addMarkovChainListener(chainListener);

        MCMCMCRunner[] threads = new MCMCMCRunner[chains.length];
        for (int i = 0; i < chains.length; i++) {
            threads[i] = new MCMCMCRunner(chains[i], mcmcmcOptions.getSwapChainsEvery(), getChainLength(), false);
            if (!DEBUG_IN_SERIES) {
                threads[i].start();
            }
        }

        while (chains[coldChain].getCurrentLength() < getChainLength()) {

            if (DEBUG_IN_SERIES) {
                for (int i = 0; i < chains.length; ++i) {
                    threads[i].runSubChain();
                }
            } else {
                // wait for all the threads to complete their alloted chain length
                boolean allDone;
                do {
                    try {
                        Thread.sleep(10);
                    } catch (InterruptedException e) {
                        //
                    }

                    allDone = true;
                    for (int i = 0; i < chains.length; i++) {
                        if (!threads[i].isChainDone()) {
                            allDone = false;
                        }
                    }
                } while (!allDone);
            }

            if (chains[coldChain].getCurrentLength() < getChainLength()) {
                int oldColdChain = coldChain;

                // attempt to swap two or more chains' temperatures
                coldChain = scheme.swapChainTemperatures(coldChain);

                // if the cold chain was involved in a swap then we need to change the
                // listener that does the logging and the destinations for the coldChainLoggers.
                if (coldChain != oldColdChain) {

                    chains[oldColdChain].removeMarkovChainListener(chainListener);

                    // Set the new cold chain's loggers with the formatters (destinations) of
                    // the original cold chain
                    for (int i = 0; i < mcLoggers[coldChain].length; i++) {
                        mcLoggers[coldChain][i].setFormatters(logFormatters[i]);
                    }

                    // Set the old cold chain to have null log formatters...
                    for (int i = 0; i < mcLoggers[oldColdChain].length; i++) {
                        mcLoggers[oldColdChain][i].setFormatters(Collections.EMPTY_LIST);
                    }

                    chains[coldChain].addMarkovChainListener(chainListener);

                }

                for (int i = 0; i < chains.length; i++) {
                    threads[i].continueChain();
                }
            }

        }

        finish();

        timer.stop();
    }

    private void runChains(long length, boolean disableCoerce) {

        Thread[] threads = new Thread[chains.length];
        for (int i = 0; i < chains.length; i++) {
            threads[i] = new MCMCMCRunner(chains[i], length, length, false);
            threads[i].start();
        }

        // wait for all threads collected to die
        for (int i = 0; i < chains.length; i++) {
            // wait doggedly for thread to die
            while (threads[i].isAlive()) {
                try {
                    threads[i].join();
                } catch (InterruptedException ie) {
                    // DO NOTHING
                }
            }

        }
    }

    private void resetChains() {

        for (MarkovChain chain : chains) {
            chain.reset();
        }
    }

    /**
     * cleans up when the chain finishes (possibly early).
     */
    private void finish() {

        NumberFormatter formatter = new NumberFormatter(8);

        MCLogger[] loggers = mcLoggers[coldChain];
        for (MCLogger logger : loggers) {
            logger.log(currentState);
            logger.stopLogging();
        }

        System.out.println();
        System.out.println("Time taken: " + timer.toString());
        System.out.println();
        System.out.println(scheme.getReport());

        if (showOperatorAnalysis) {
            OperatorAnalysisPrinter.showOperatorAnalysis(System.out, schedules[coldChain], mcmcOptions.useAdaptation(), false);
            System.out.println();
        }
    }

    private final MarkovChainListener chainListener = new MarkovChainListener() {

        // MarkovChainListener interface *******************************************
        // for receiving messages from subordinate MarkovChain

        /**
         * Called to update the current model keepEvery states.
         */
        public synchronized void currentState(long state, MarkovChain markovChain, Model currentModel) {

            currentState = state;

            if (state % 1000 == 0) {
                NumberFormatter formatter = new NumberFormatter(8);
                formatter.setPadding(false);

                //System.out.print("State " + currentState + ": ");
//                if (coldChain != 0) {
//                    MCMCCriterion acceptor = ((MCMCCriterion) chains[coldChain].getAcceptor());
//                    int rank = acceptor.getRank();
//                    System.err.println("here");
//                }
                for (int i = 0; i < chains.length; i++) {
                    String score;
                    if (i == coldChain) {
                        score = "[" + formatter.format(chains[i].getCurrentScore()) + "]";
                    } else {
                        score = formatter.format(chains[i].getCurrentScore());
                    }
                    score += " ";
                    System.out.print(formatter.formatToFieldWidth(score, 12));
                }
                System.out.println();
            }

            MCLogger[] loggers = mcLoggers[coldChain];
            for (MCLogger logger : loggers) {
                logger.log(state);
            }

        }

        /**
         * Called when a new best posterior state is found.
         */
        public synchronized void bestState(long state, MarkovChain markovChain, Model bestModel) {
            currentState = state;
        }

        /**
         * cleans up when the chain finishes (possibly early).
         */
        public synchronized void finished(long chainLength, MarkovChain markovChain) {
        }

    };

    public int getColdChain() {
        return coldChain;
    }
    
    /**
     * @return the likelihood function.
     */
    public Likelihood getLikelihood() {
        return chains[coldChain].getLikelihood();
    }

    /**
     * @return the timer.
     */
    public dr.util.Timer getTimer() {
        return timer;
    }

    /**
     * @return the length of this analysis.
     */
    public final long getChainLength() {
        return mcmcOptions.getChainLength();
    }

    // TRANSIENT PUBLIC METHODS *****************************************

    /**
     * @return the current state of the MCMC analysis.
     */
    public final long getCurrentState() {
        return currentState;
    }

    /**
     * @return the progress (0 to 1) of the MCMC analysis.
     */
    public final double getProgress() {
        return (double) currentState / (double) mcmcOptions.getChainLength();
    }

    /**
     * Requests that the MCMC chain stop prematurely.
     */
    public void pleaseStop() {
        for (MarkovChain chain : chains) {
            chain.pleaseStop();
        }
    }

    public void setShowOperatorAnalysis(boolean soa) {
        showOperatorAnalysis = soa;
    }

    // PRIVATE TRANSIENTS

    private final MCMCOptions mcmcOptions;
    private final MCMCMCOptions mcmcmcOptions;

    private boolean showOperatorAnalysis = true;
    private final dr.util.Timer timer = new dr.util.Timer();
    private long currentState = 0;

    private final MarkovChain[] chains;
    private final MCLogger[][] mcLoggers;
    private final OperatorSchedule[] schedules;
    private int coldChain;

    private final ParallelTempering scheme;
}

