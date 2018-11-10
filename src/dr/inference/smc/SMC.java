/*
 * MCMC.java
 *
 * Copyright (c) 2002-2015 Alexei Drummond, Andrew Rambaut and Marc Suchard
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
 */

package dr.inference.smc;

import dr.app.checkpoint.BeastCheckpointer;
import dr.inference.loggers.LogColumn;
import dr.inference.loggers.Loggable;
import dr.inference.loggers.Logger;
import dr.inference.markovchain.MarkovChain;
import dr.inference.markovchain.MarkovChainListener;
import dr.inference.mcmc.MCMCCriterion;
import dr.inference.mcmc.MCMCOptions;
import dr.inference.model.Likelihood;
import dr.inference.model.Model;
import dr.inference.operators.*;
import dr.inference.state.Factory;
import dr.inference.state.StateLoader;
import dr.inference.state.StateSaver;
import dr.util.Identifiable;
import dr.util.NumberFormatter;
import dr.xml.Spawnable;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;

/**
 * A class that runs short MCMC chains for each of a set of particles as
 * part of a sequential Monte Carlo (SMC) sampler.
 *
 * @author Andrew Rambaut
 * @version $Id:$
 */
public class SMC implements Identifiable, Spawnable, Loggable {

    public SMC(String id, List<BeastCheckpointer> particleStates) {
        this.id = id;
        this.particleStates.addAll(particleStates);
    }

    /**
     * Must be called before calling chain.
     *
     * @param options    the options for this SMC analysis
     * @param schedule   operator schedule to be used in chain.
     * @param likelihood the likelihood for this SMC
     * @param loggers    an array of loggers to record output of this SMC run
     */
    public void init(
            SMCOptions options,
            Likelihood likelihood,
            OperatorSchedule schedule,
            Logger[] loggers) {

        MCMCCriterion criterion = new MCMCCriterion();

        // full evaluation tests and operator adaptation are off as these are multiple short runs.
        // Operator tuning will have already been done.
        mc = new MarkovChain(likelihood, schedule, criterion, 0, 0, 0, false);

        this.options = options;
        this.loggers = loggers;
        this.schedule = schedule;

        //initialize transients
        currentState = 0;

        // States are saved at the end of each particle's run
//        if (Factory.INSTANCE != null) {
//            for (MarkovChainListener listener : Factory.INSTANCE.getStateSaverChainListeners()) {
//                mc.addMarkovChainListener(listener);
//            }
//        }

    }

    public MarkovChain getMarkovChain() {
        return mc;
    }

    public Logger[] getLoggers() {
        return loggers;
    }

    public OperatorSchedule getOperatorSchedule() {
        return schedule;
    }

    public void run() {
        chain();
    }

    /**
     * This method actually initiates the MCMC analysis.
     */
    public void chain() {

        currentState = 0;

        timer.start();

        if (loggers != null) {
            for (Logger logger : loggers) {
                logger.startLogging();
            }
        }

        mc.addMarkovChainListener(chainListener);

        for (BeastCheckpointer particleState : particleStates) {
            // Don't need the savedLnL - it won't be there
            particleState.loadState(mc, new double[1]);

            // reset the current chain length to 0
            mc.setCurrentLength(0);

            mc.runChain(options.getChainLength(), true);

            // Save state to file...
            particleState.saveState(mc, mc.getCurrentLength(), mc.getCurrentScore());
        }

        mc.terminateChain();

        mc.removeMarkovChainListener(chainListener);

        timer.stop();
    }

    @Override
    public LogColumn[] getColumns() {
        return new LogColumn[] { new LogColumn() {
            @Override
            public void setLabel(String label) {
            }

            @Override
            public String getLabel() {
                return "time";
            }

            @Override
            public void setMinimumWidth(int minimumWidth) {

            }

            @Override
            public int getMinimumWidth() {
                return 0;
            }

            @Override
            public String getFormatted() {
                return Double.toString(getTimer().toSeconds());
            }
        }   };
    }

    private final MarkovChainListener chainListener = new MarkovChainListener() {

        // MarkovChainListener interface *******************************************
        // for receiving messages from subordinate MarkovChain

        /**
         * Called to update the current model keepEvery states.
         */
        @Override
        public void currentState(long state, MarkovChain markovChain, Model currentModel) {

            currentState = state;

            if (loggers != null) {
                for (Logger logger : loggers) {
                    logger.log(state);
                }
            }
        }

        /**
         * Called when a new new best posterior state is found.
         */
        @Override
        public void bestState(long state, MarkovChain markovChain, Model bestModel) { }

        /**
         * cleans up when the chain finishes (possibly early).
         */
        @Override
        public void finished(long chainLength, MarkovChain markovChain) {
            currentState = chainLength;

            if (loggers != null) {
                for (Logger logger : loggers) {
                    logger.log(currentState);
                    logger.stopLogging();
                }
            }
        }

    };

    /**
     * @return the likelihood function.
     */
    public Likelihood getLikelihood() {
        return mc.getLikelihood();
    }

    /**
     * @return the timer.
     */
    public dr.util.Timer getTimer() {
        return timer;
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
        return (double) currentState / (double) options.getChainLength();
    }

    public boolean getSpawnable() {
        return true;
    }

    private boolean spawnable = true;

    //PRIVATE METHODS *****************************************

    public String getId() {
        return id;
    }

    public void setId(String id) {
        this.id = id;
    }

    // PRIVATE TRANSIENTS

    private final dr.util.Timer timer = new dr.util.Timer();
    private long currentState = 0;
    private final NumberFormatter formatter = new NumberFormatter(8);

    /**
     * this markov chain does most of the work.
     */
    private MarkovChain mc;

    /**
     * the options of this MCMC analysis
     */
    private SMCOptions options;

    private final List<BeastCheckpointer> particleStates = new ArrayList<BeastCheckpointer>();
    private Logger[] loggers;
    private OperatorSchedule schedule;

    private String id = null;
}
