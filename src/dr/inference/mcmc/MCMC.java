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

package dr.inference.mcmc;

import dr.inference.state.Factory;
import dr.inference.state.StateLoader;
import dr.inference.loggers.LogColumn;
import dr.inference.loggers.Loggable;
import dr.inference.loggers.Logger;
import dr.inference.markovchain.MarkovChain;
import dr.inference.markovchain.MarkovChainListener;
import dr.inference.model.Likelihood;
import dr.inference.model.Model;
import dr.inference.operators.*;
import dr.util.Identifiable;
import dr.util.NumberFormatter;
import dr.xml.Spawnable;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.PrintStream;

/**
 * An MCMC analysis that estimates parameters of a probabilistic model.
 *
 * @author Alexei Drummond
 * @author Andrew Rambaut
 * @version $Id: MCMC.java,v 1.41 2005/07/11 14:06:25 rambaut Exp $
 */
public class MCMC implements Identifiable, Spawnable, Loggable {

    public MCMC(String id) {
        this.id = id;
    }

    /**
     * Must be called before calling chain.
     *
     * @param options    the options for this MCMC analysis
     * @param schedule   operator schedule to be used in chain.
     * @param likelihood the likelihood for this MCMC
     * @param loggers    an array of loggers to record output of this MCMC run
     */
    public void init(
            MCMCOptions options,
            Likelihood likelihood,
            OperatorSchedule schedule,
            Logger[] loggers) {

        MCMCCriterion criterion = new MCMCCriterion();
        criterion.setTemperature(options.getTemperature());

        mc = new MarkovChain(likelihood, schedule, criterion,
                options.getFullEvaluationCount(), options.minOperatorCountForFullEvaluation(),
                options.getEvaluationTestThreshold(),
                options.useCoercion());

        this.options = options;
        this.loggers = loggers;
        this.schedule = schedule;

        //initialize transients
        currentState = 0;

        if (Factory.INSTANCE != null) {
            for (MarkovChainListener listener : Factory.INSTANCE.getStateSaverChainListeners()) {
                mc.addMarkovChainListener(listener);
            }
        }

    }

    /**
     * Must be called before calling chain.
     *
     * @param chainlength chain length
     * @param likelihood the likelihood for this MCMC
     * @param operators  an array of MCMC operators
     * @param loggers    an array of loggers to record output of this MCMC run
     */
    public void init(long chainlength,
                     Likelihood likelihood,
                     MCMCOperator[] operators,
                     Logger[] loggers) {

        MCMCOptions options = new MCMCOptions(chainlength);
        MCMCCriterion criterion = new MCMCCriterion();
        criterion.setTemperature(1);
        OperatorSchedule schedule = new SimpleOperatorSchedule();
        for (MCMCOperator operator : operators) schedule.addOperator(operator);

        init(options, likelihood, schedule, loggers);
    }

    public MarkovChain getMarkovChain() {
        return mc;
    }

    public Logger[] getLoggers() {
        return loggers;
    }

    public MCMCOptions getOptions() {
        return options;
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

        stopping = false;
        currentState = 0;

        timer.start();

        if (loggers != null) {
            for (Logger logger : loggers) {
                logger.startLogging();
            }
        }

        if (!stopping) {

            long loadedState = 0;

            if (Factory.INSTANCE != null) {
                StateLoader initialStateLoader = Factory.INSTANCE.getInitialStateLoader();
                if (initialStateLoader != null) {
                    double[] savedLnL = new double[1];

                    initialStateLoader.loadState(mc, savedLnL);

                    mc.setCurrentLength(loadedState);

                    double lnL = mc.evaluate();

                    initialStateLoader.checkLoadState(savedLnL[0], lnL);
                }
            }

            mc.addMarkovChainListener(chainListener);

            long chainLength = getChainLength();

            //this also potentially gets the new coercionDelay of a possibly increased chain length
            final long coercionDelay = getCoercionDelay();

            //assume that dumped state has passed the coercionDelay
            //TODO: discuss whether we want to dump the coercionDelay or chainLength to file
            if (coercionDelay > loadedState) {
                mc.runChain(coercionDelay - loadedState, true);
                chainLength -= coercionDelay;
            }

            //if (coercionDelay > 0) {
            // Run the chain for coercionDelay steps with coercion disabled
            //mc.runChain(coercionDelay, true);
            //chainLength -= coercionDelay;

            // reset operator acceptance levels
            //GB: we are now restoring these; commenting out for now
                /*for (int i = 0; i < schedule.getOperatorCount(); i++) {
                    schedule.getOperator(i).reset();
                }*/
            //}

            mc.runChain(chainLength, false);

            mc.terminateChain();

            mc.removeMarkovChainListener(chainListener);

        }
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
            // OperatorAnalysisPrinter class can do the job now
            if (showOperatorAnalysis) {
                OperatorAnalysisPrinter.showOperatorAnalysis(System.out, getOperatorSchedule(), options.useCoercion());
            }

            if (operatorAnalysisFile != null) {
                try {
                    PrintStream out = new PrintStream(new FileOutputStream(operatorAnalysisFile));
                    OperatorAnalysisPrinter.showOperatorAnalysis(out, getOperatorSchedule(), options.useCoercion());
                    out.flush();
                    out.close();
                } catch (IOException e) {
                    e.printStackTrace();
                }
            }

            // How should premature finish be flagged?
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

    /**
     * @return the length of this analysis.
     */
    public final long getChainLength() {
        return options.getChainLength();
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

    /**
     * @return true if this MCMC is currently adapting the operators.
     */
    public final boolean isAdapting() {
        return isAdapting;
    }

    /**
     * Requests that the MCMC chain stop prematurely.
     */
    public void pleaseStop() {
        stopping = true;
        mc.pleaseStop();
    }

    /**
     * @return true if Markov chain is stopped
     */
    public boolean isStopped() {
        return mc.isStopped();
    }

    public boolean getSpawnable() {
        return spawnable;
    }

    private boolean spawnable = true;

    public void setSpawnable(boolean spawnable) {
        this.spawnable = spawnable;
    }


    //PRIVATE METHODS *****************************************
    protected long getCoercionDelay() {

        long delay = options.getCoercionDelay();
        if (delay < 0) {
            delay = (long)(options.getChainLength() / 100);
        }
        if (options.useCoercion()) return delay;

        for (int i = 0; i < schedule.getOperatorCount(); i++) {
            MCMCOperator op = schedule.getOperator(i);

            if (op instanceof CoercableMCMCOperator) {
                if (((CoercableMCMCOperator) op).getMode() == CoercionMode.COERCION_ON) return delay;
            }
        }

        return -1;
    }

    public void setShowOperatorAnalysis(boolean soa) {
        showOperatorAnalysis = soa;
    }


    public void setOperatorAnalysisFile(File operatorAnalysisFile) {
        this.operatorAnalysisFile = operatorAnalysisFile;
    }

    public String getId() {
        return id;
    }

    public void setId(String id) {
        this.id = id;
    }

    // PRIVATE TRANSIENTS

    //private FileLogger operatorLogger = null;
    protected final boolean isAdapting = true;
    protected boolean stopping = false;
    protected boolean showOperatorAnalysis = true;
    protected File operatorAnalysisFile = null;
    protected final dr.util.Timer timer = new dr.util.Timer();
    protected long currentState = 0;
    //private int stepsPerReport = 1000;
    protected final NumberFormatter formatter = new NumberFormatter(8);

    /**
     * this markov chain does most of the work.
     */
    protected MarkovChain mc;

    /**
     * the options of this MCMC analysis
     */
    protected MCMCOptions options;

    protected Logger[] loggers;
    protected OperatorSchedule schedule;

    private String id = null;
}
