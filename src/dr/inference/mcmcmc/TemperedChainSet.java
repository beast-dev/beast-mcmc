/*
 * TemperedChainSet.java
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
import dr.inference.mcmc.ChainContext;
import dr.inference.mcmc.MCMC;
import dr.inference.mcmc.MCMCCriterion;
import dr.inference.operators.OperatorSchedule;

import java.io.File;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.concurrent.Callable;
import java.util.concurrent.Semaphore;

/**
 * Metropolis-coupled MCMC (parallel tempering): a cold chain and one or more heated chains
 * that run in parallel between swap attempts. At each swap the chains' temperatures (and
 * operator tunings) are exchanged, and the log destinations follow the cold temperature so
 * the output is always written by whichever chain is currently cold.
 *
 * @author Andrew Rambaut
 */
public class TemperedChainSet extends AbstractChainSet {

    public TemperedChainSet(int index, MCMC[] mcmcs, ChainContext[] contexts, ChainContext driverContext,
                            MCMCMCOptions options) {
        super(index, mcmcs, contexts, driverContext);

        this.options = options;

        double[] temperatures = options.getChainTemperatures();
        if (temperatures[0] != 1.0) {
            throw new IllegalArgumentException("The first chain in the array should be cold (temperature = 1.0)");
        }
        if (temperatures.length != mcmcs.length) {
            throw new IllegalArgumentException("There should be one temperature per chain");
        }

        chains = new MarkovChain[mcmcs.length];
        schedules = new OperatorSchedule[mcmcs.length];
        mcLoggers = new MCLogger[mcmcs.length][];

        for (int i = 0; i < mcmcs.length; i++) {
            chains[i] = mcmcs[i].getMarkovChain();
            schedules[i] = mcmcs[i].getOperatorSchedule();

            MCMCCriterion acceptor = (MCMCCriterion) chains[i].getAcceptor();
            acceptor.setTemperature(temperatures[i]);
            acceptor.setRank(i);

            // Every chain logs, but only the cold chain's loggers have destinations. The
            // destinations are moved to whichever chain becomes cold at a swap.
            Logger[] loggers = mcmcs[i].getLoggers();
            mcLoggers[i] = new MCLogger[loggers == null ? 0 : loggers.length];
            for (int j = 0; j < mcLoggers[i].length; j++) {
                mcLoggers[i][j] = (MCLogger) loggers[j];
            }
        }

        // only the first chain was built with real output files
        operatorAnalysisFile = mcmcs[0].getOperatorAnalysisFile();
        operatorAnalysisCSVFile = mcmcs[0].isOperatorAnalysisCSVFile();

        coldChain = 0;
        scheme = options.getSwapScheme().factory(chains, schedules, options);
    }

    @Override
    public void start(Semaphore permits) {
        super.start(permits);

        coldFormatters = new ArrayList<>();
        for (MCLogger logger : mcLoggers[coldChain]) {
            coldFormatters.add(logger.getFormatters());
        }
        for (int i = 0; i < mcmcs.length; i++) {
            if (i != coldChain) {
                for (MCLogger logger : mcLoggers[i]) {
                    logger.setFormatters(Collections.<LogFormatter>emptyList());
                }
            }
        }

        List<Callable<Boolean>> tasks = new ArrayList<>();
        for (final MCMC mcmc : mcmcs) {
            tasks.add(() -> {
                mcmc.start();
                return true;
            });
        }
        runOnChains(tasks);
    }

    @Override
    public boolean step() {
        final long swapEvery = options.getSwapChainsEvery();

        List<Callable<Boolean>> tasks = new ArrayList<>();
        for (final MCMC mcmc : mcmcs) {
            tasks.add(() -> mcmc.runSegment(swapEvery));
        }

        // runOnChains waits for every chain, so the swap below sees their final state
        boolean more = true;
        for (boolean chainHasMore : runOnChains(tasks)) {
            more &= chainHasMore;
        }

        if (more) {
            int oldColdChain = coldChain;
            coldChain = scheme.swapChainTemperatures(coldChain);
            if (coldChain != oldColdChain) {
                moveLogDestinations(oldColdChain, coldChain);
            }
        }

        return more;
    }

    private void moveLogDestinations(int from, int to) {
        for (int j = 0; j < mcLoggers[to].length; j++) {
            mcLoggers[to][j].setFormatters(coldFormatters.get(j));
        }
        for (MCLogger logger : mcLoggers[from]) {
            logger.setFormatters(Collections.<LogFormatter>emptyList());
        }
    }

    @Override
    public void finish() {
        // The operator tunings and counts move with the temperatures, so the cold chain's
        // schedule holds the cold statistics. Only report those.
        for (int i = 0; i < mcmcs.length; i++) {
            if (i != coldChain) {
                mcmcs[i].setShowOperatorAnalysis(false);
                mcmcs[i].setOperatorAnalysisFile(null, false);
            } else {
                mcmcs[i].setOperatorAnalysisFile(operatorAnalysisFile, operatorAnalysisCSVFile);
            }
        }

        List<Callable<Boolean>> tasks = new ArrayList<>();
        for (final MCMC mcmc : mcmcs) {
            tasks.add(() -> {
                mcmc.finish();
                return true;
            });
        }
        runOnChains(tasks);

        System.out.println();
        System.out.println(getReport());
    }

    @Override
    public MCMC getColdMCMC() {
        return mcmcs[coldChain];
    }

    @Override
    public MarkovChain getColdChain() {
        return chains[coldChain];
    }

    @Override
    public long getCurrentState() {
        return chains[coldChain].getCurrentLength();
    }

    @Override
    public String getStatus() {
        return String.format("swap acceptance %.3f", scheme.getStatistics().getAcceptanceRate());
    }

    public String getReport() {
        return scheme.getReport();
    }

    private final MCMCMCOptions options;
    private final MarkovChain[] chains;
    private final OperatorSchedule[] schedules;
    private final MCLogger[][] mcLoggers;
    private final ParallelTempering scheme;

    private final File operatorAnalysisFile;
    private final boolean operatorAnalysisCSVFile;

    private List<List<LogFormatter>> coldFormatters;
    private volatile int coldChain;
}
