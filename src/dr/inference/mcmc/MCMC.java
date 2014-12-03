/*
 * MCMC.java
 *
 * Copyright (C) 2002-2009 Alexei Drummond and Andrew Rambaut
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
 * BEAST is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with BEAST; if not, write to the
 * Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
 * Boston, MA  02110-1301  USA
 */

package dr.inference.mcmc;

import dr.inference.loggers.Logger;
import dr.inference.markovchain.MarkovChain;
import dr.inference.markovchain.MarkovChainDelegate;
import dr.inference.markovchain.MarkovChainListener;
import dr.inference.model.Likelihood;
import dr.inference.model.Model;
import dr.inference.operators.*;
import dr.inference.prior.Prior;
import dr.util.Identifiable;
import dr.util.NumberFormatter;
import dr.xml.Spawnable;

import java.io.*;
import java.util.Calendar;

/**
 * An MCMC analysis that estimates parameters of a probabilistic model.
 *
 * @author Alexei Drummond
 * @author Andrew Rambaut
 * @version $Id: MCMC.java,v 1.41 2005/07/11 14:06:25 rambaut Exp $
 */
public class MCMC implements Identifiable, Spawnable {

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

        init(options, likelihood, Prior.UNIFORM_PRIOR, schedule, loggers, new MarkovChainDelegate[0]);
    }

    /**
     * Must be called before calling chain.
     *
     * @param options    the options for this MCMC analysis
     * @param schedule   operator schedule to be used in chain.
     * @param likelihood the likelihood for this MCMC
     * @param loggers    an array of loggers to record output of this MCMC run
     * @param delegates    an array of delegates to handle tasks related to the MCMC
     */
    public void init(
            MCMCOptions options,
            Likelihood likelihood,
            OperatorSchedule schedule,
            Logger[] loggers,
            MarkovChainDelegate[] delegates) {

        init(options, likelihood, Prior.UNIFORM_PRIOR, schedule, loggers, delegates);
    }

    /**
     * Must be called before calling chain.
     *
     * @param options    the options for this MCMC analysis
     * @param prior      the prior disitrbution on the model parameters.
     * @param schedule   operator schedule to be used in chain.
     * @param likelihood the likelihood for this MCMC
     * @param loggers    an array of loggers to record output of this MCMC run
     */
    public void init(
            MCMCOptions options,
            Likelihood likelihood,
            Prior prior,
            OperatorSchedule schedule,
            Logger[] loggers) {

        init(options, likelihood, prior, schedule, loggers, new MarkovChainDelegate[0]);
    }

    /**
     * Must be called before calling chain.
     *
     * @param options    the options for this MCMC analysis
     * @param prior      the prior disitrbution on the model parameters.
     * @param schedule   operator schedule to be used in chain.
     * @param likelihood the likelihood for this MCMC
     * @param loggers    an array of loggers to record output of this MCMC run
     * @param delegates    an array of delegates to handle tasks related to the MCMC
     */
    public void init(
            MCMCOptions options,
            Likelihood likelihood,
            Prior prior,
            OperatorSchedule schedule,
            Logger[] loggers,
            MarkovChainDelegate[] delegates) {

        MCMCCriterion criterion = new MCMCCriterion();
        criterion.setTemperature(options.getTemperature());

        mc = new MarkovChain(prior, likelihood, schedule, criterion,
                options.getFullEvaluationCount(), options.minOperatorCountForFullEvaluation(),
                options.getEvaluationTestThreshold(),
                options.useCoercion());

        this.options = options;
        this.loggers = loggers;
        this.schedule = schedule;

        //initialize transients
        currentState = 0;

        // Does not seem to be in use (JH)
/*
        stepsPerReport = 1;
        while ((getChainLength() / stepsPerReport) > 1000) {
            stepsPerReport *= 2;
        }*/

        for(MarkovChainDelegate delegate : delegates) {
            delegate.setup(options, schedule, mc);
        }
        this.delegates = delegates;
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

        init(options, likelihood, Prior.UNIFORM_PRIOR, schedule, loggers);
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

    // Experimental
    public final static boolean TEST_CLONING = false;

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
            mc.addMarkovChainListener(chainListener);

            for(MarkovChainDelegate delegate : delegates) {
                mc.addMarkovChainDelegate(delegate);
            }

            long chainLength = getChainLength();

            final long coercionDelay = getCoercionDelay();

            if (coercionDelay > 0) {
                // Run the chain for coercionDelay steps with coercion disabled
                mc.runChain(coercionDelay, true);
                chainLength -= coercionDelay;

                // reset operator acceptance levels
                for (int i = 0; i < schedule.getOperatorCount(); i++) {
                    schedule.getOperator(i).reset();
                }
            }

            if (TEST_CLONING) {
                // TEST Code for cloning the MarkovChain to a file to distribute amongst processors
                for(MarkovChainDelegate delegate : delegates) {
                    mc.removeMarkovChainDelegate(delegate);
                }
                mc.removeMarkovChainListener(chainListener);

                // Write the MarkovChain out and back in again...
                writeMarkovChainToFile(new File("beast.clone"), mc);
                mc = readMarkovChainFromFile(new File("beast.clone"));

                mc.addMarkovChainListener(chainListener);

                for(MarkovChainDelegate delegate : delegates) {
                    mc.addMarkovChainDelegate(delegate);
                }
                // TEST Code end
            }

            mc.runChain(chainLength, false);

            mc.terminateChain();

            mc.removeMarkovChainListener(chainListener);

            for(MarkovChainDelegate delegate : delegates) {
                mc.removeMarkovChainDelegate(delegate);
            }
        }
        timer.stop();
    }

    protected boolean writeMarkovChainToFile(File file, MarkovChain mc) {
        OutputStream fileOut = null;
        try {
            fileOut = new FileOutputStream(file);
            ObjectOutputStream out = new ObjectOutputStream(fileOut);
            out.writeObject(mc);
            out.close();
            fileOut.close();
        } catch (IOException ioe) {
            System.err.println("Unable to write file: " + ioe.getMessage());
        }
        return true;
    }

    protected MarkovChain readMarkovChainFromFile(File file) {
        MarkovChain mc = null;
        try {
            FileInputStream fileIn =
                    new FileInputStream(file);
            ObjectInputStream in = new ObjectInputStream(fileIn);
            mc = (MarkovChain) in.readObject();
            in.close();
            fileIn.close();
        } catch (IOException ioe) {
            System.err.println("Unable to read file: " + ioe.getMessage());
        } catch (ClassNotFoundException cnfe) {
            System.err.println("Unable to read file: " + cnfe.getMessage());
            cnfe.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
        }

        return mc;
    }

    protected final MarkovChainListener chainListener = new MarkovChainListener() {

        // MarkovChainListener interface *******************************************
        // for receiving messages from subordinate MarkovChain

        /**
         * Called to update the current model keepEvery states.
         */
        public void currentState(long state, Model currentModel) {

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
        public void bestState(long state, Model bestModel) {
            currentState = state;
        }

        /**
         * cleans up when the chain finishes (possibly early).
         */
        public void finished(long chainLength) {
            currentState = chainLength;

            if (loggers != null) {
                for (Logger logger : loggers) {
                    logger.log(currentState);
                    logger.stopLogging();
                }
            }
            // OperatorAnalysisPrinter class can do the job now
            if (showOperatorAnalysis) {
                showOperatorAnalysis(System.out);
            }

            if (operatorAnalysisFile != null) {
                try {
                    PrintStream out = new PrintStream(new FileOutputStream(operatorAnalysisFile));
                    showOperatorAnalysis(out);
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
     * Writes ano operator analysis to the provided print stream
     *
     * @param out the print stream to write operator analysis to
     */
    private void showOperatorAnalysis(PrintStream out) {
        out.println();
        out.println("Operator analysis");
        out.println(formatter.formatToFieldWidth("Operator", 50) +
                formatter.formatToFieldWidth("Tuning", 9) +
                formatter.formatToFieldWidth("Count", 11) +
                formatter.formatToFieldWidth("Time", 9) +
                formatter.formatToFieldWidth("Time/Op", 9) +
                formatter.formatToFieldWidth("Pr(accept)", 11) +
                (options.useCoercion() ? "" : " Performance suggestion"));

        for (int i = 0; i < schedule.getOperatorCount(); i++) {

            final MCMCOperator op = schedule.getOperator(i);
            if (op instanceof JointOperator) {
                JointOperator jointOp = (JointOperator) op;
                for (int k = 0; k < jointOp.getNumberOfSubOperators(); k++) {
                    out.println(formattedOperatorName(jointOp.getSubOperatorName(k))
                            + formattedParameterString(jointOp.getSubOperator(k))
                            + formattedCountString(op)
                            + formattedTimeString(op)
                            + formattedTimePerOpString(op)
                            + formattedProbString(jointOp)
                            + (options.useCoercion() ? "" : formattedDiagnostics(jointOp, MCMCOperator.Utils.getAcceptanceProbability(jointOp)))
                    );
                }
            } else {
                out.println(formattedOperatorName(op.getOperatorName())
                        + formattedParameterString(op)
                        + formattedCountString(op)
                        + formattedTimeString(op)
                        + formattedTimePerOpString(op)
                        + formattedProbString(op)
                        + (options.useCoercion() ? "" : formattedDiagnostics(op, MCMCOperator.Utils.getAcceptanceProbability(op)))
                );
            }

        }
        out.println();
    }

    private String formattedOperatorName(String operatorName) {
        return formatter.formatToFieldWidth(operatorName, 50);
    }

    private String formattedParameterString(MCMCOperator op) {
        String pString = "        ";
        if (op instanceof CoercableMCMCOperator && ((CoercableMCMCOperator) op).getMode() != CoercionMode.COERCION_OFF) {
            pString = formatter.formatToFieldWidth(formatter.formatDecimal(((CoercableMCMCOperator) op).getRawParameter(), 3), 8);
        }
        return pString;
    }

    private String formattedCountString(MCMCOperator op) {
        final int count = op.getCount();
        return formatter.formatToFieldWidth(Integer.toString(count), 10) + " ";
    }

    private String formattedTimeString(MCMCOperator op) {
        final long time = op.getTotalEvaluationTime();
        return formatter.formatToFieldWidth(Long.toString(time), 8) + " ";
    }

    private String formattedTimePerOpString(MCMCOperator op) {
        final double time = op.getMeanEvaluationTime();
        return formatter.formatToFieldWidth(formatter.formatDecimal(time, 2), 8) + " ";
    }

    private String formattedProbString(MCMCOperator op) {
        final double acceptanceProb = MCMCOperator.Utils.getAcceptanceProbability(op);
        return formatter.formatToFieldWidth(formatter.formatDecimal(acceptanceProb, 4), 11) + " ";
    }

    private String formattedDiagnostics(MCMCOperator op, double acceptanceProb) {

        String message = "good";
        if (acceptanceProb < op.getMinimumGoodAcceptanceLevel()) {
            if (acceptanceProb < (op.getMinimumAcceptanceLevel() / 10.0)) {
                message = "very low";
            } else if (acceptanceProb < op.getMinimumAcceptanceLevel()) {
                message = "low";
            } else message = "slightly low";

        } else if (acceptanceProb > op.getMaximumGoodAcceptanceLevel()) {
            double reallyHigh = 1.0 - ((1.0 - op.getMaximumAcceptanceLevel()) / 10.0);
            if (acceptanceProb > reallyHigh) {
                message = "very high";
            } else if (acceptanceProb > op.getMaximumAcceptanceLevel()) {
                message = "high";
            } else message = "slightly high";
        }

        String performacsMsg;
        if (op instanceof GibbsOperator) {
            performacsMsg = "none (Gibbs operator)";
        } else {
            final String suggestion = op.getPerformanceSuggestion();
            performacsMsg = message + "\t" + suggestion;
        }

        return performacsMsg;
    }

    /**
     * @return the prior of this MCMC analysis.
     */
    public Prior getPrior() {
        return mc.getPrior();
    }

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
    private MarkovChainDelegate[] delegates;

    private String id = null;
}

