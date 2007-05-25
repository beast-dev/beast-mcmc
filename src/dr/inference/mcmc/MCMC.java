/*
 * MCMC.java
 *
 * Copyright (C) 2002-2006 Alexei Drummond and Andrew Rambaut
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

import dr.inference.loggers.Logger;
import dr.inference.markovchain.MarkovChain;
import dr.inference.markovchain.MarkovChainListener;
import dr.inference.model.Likelihood;
import dr.inference.model.Model;
import dr.inference.model.CompoundLikelihood;
import dr.inference.operators.CoercableMCMCOperator;
import dr.inference.operators.MCMCOperator;
import dr.inference.operators.OperatorSchedule;
import dr.inference.prior.Prior;
import dr.util.Identifiable;
import dr.util.NumberFormatter;
import dr.xml.*;

import java.util.ArrayList;

/**
 * An MCMC analysis that estimates parameters of a probabilistic model.
 *
 * @author Alexei Drummond
 * @author Andrew Rambaut
 *
 * @version $Id: MCMC.java,v 1.41 2005/07/11 14:06:25 rambaut Exp $
 */
public class MCMC implements Runnable, Identifiable {

    public MCMC(String id) {
        this.id = id;
    }

    /**
     * Must be called before calling chain.
     * @param options the options for this MCMC analysis
     * @param prior the prior disitrbution on the model parameters.
     * @param schedule operator schedule to be used in chain.
     */
    public void init(
            MCMCOptions options,
            Likelihood likelihood,
            Prior prior,
            OperatorSchedule schedule,
            Logger[] loggers) {

        MCMCCriterion criterion = new MCMCCriterion();
        criterion.setTemperature(options.getTemperature());

        mc = new MarkovChain(prior, likelihood, schedule, criterion, options.useCoercion());

        this.options = options;
        this.loggers = loggers;
        this.schedule = schedule;

        //initialize transients
        currentState = 0;

        stepsPerReport = 1;
        while ((getChainLength() / stepsPerReport) > 1000) {
            stepsPerReport *= 2;
        }
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
     * This method actually intiates the MCMC analysis.
     */
    public void chain() {

        stopping = false;
        currentState = 0;

        timer.start();

        if (isPreBurninNeeded()) {
            int preBurnin = options.getPreBurnin();
            if (preBurnin > 0) {
                MarkovChainListener burninListener = new BurninListener(preBurnin);
                mc.addMarkovChainListener(burninListener);
                mc.chain(preBurnin, true);
                mc.removeMarkovChainListener(burninListener);
                mc.reset();
            }
        }

        if (loggers != null) {
            for (int i =0; i < loggers.length; i++) {
                loggers[i].startLogging();
            }
        }

        if (!stopping) {
            mc.addMarkovChainListener(chainListener);
            mc.chain(getChainLength(), false);
            mc.removeMarkovChainListener(chainListener);
        }
        timer.stop();
    }

    public class BurninListener implements MarkovChainListener {

        public BurninListener(int stateCount) {
            this.stateCount = stateCount;
            step = 0;
            stepSize = (double)stateCount / 60.0;
        }

        /**
         * Called to update the current model keepEvery states.
         */
        public void currentState(int state, Model currentModel) {

            if (state == 0) {
                System.out.println();
                System.out.println("Pre-burnin (" + stateCount + " states)");
                System.out.println("0              25             50             75            100");
                System.out.println("|--------------|--------------|--------------|--------------|");
                System.out.print(  "*");
                step = 1;
            }

            if (state >= (int)Math.round(step*stepSize) && step <= 60) {
                System.out.print("*");
                System.out.flush();
                step += 1;
            }
        }

        /** Called when a new new best posterior state is found. */
        public void bestState(int state, Model bestModel) {}

        /** cleans up when the chain finishes (possibly early). */
        public void finished(int chainLength) {
            System.out.println("*");
            System.out.println();
        }

        int stateCount = 0;
        double stepSize;
        int step = 0;
    }

    private MarkovChainListener chainListener = new MarkovChainListener() {

        // MarkovChainListener interface *******************************************
        // for receiving messages from subordinate MarkovChain

        /**
         * Called to update the current model keepEvery states.
         */
        public void currentState(int state, Model currentModel) {

            currentState = state;

            if (loggers != null) {
                for (int i =0; i < loggers.length; i++) {
                    loggers[i].log(state);
                }
            }
        }

        /** Called when a new new best posterior state is found. */
        public void bestState(int state, Model bestModel) {
            currentState = state;
        }

        /** cleans up when the chain finishes (possibly early). */
        public void finished(int chainLength) {
            currentState = chainLength;

            if (loggers != null) {
                for (int i =0; i < loggers.length; i++) {
                    loggers[i].log(currentState);
                    loggers[i].stopLogging();
                }
            }

            if (showOperatorAnalysis) {
                System.out.println();
                System.out.println("Operator analysis");
                System.out.println(
                        formatter.formatToFieldWidth("Operator", 30) +
                                formatter.formatToFieldWidth("", 8) +
                                formatter.formatToFieldWidth("Pr(accept)", 11) +
                                " Performance suggestion");
                for (int i =0; i < schedule.getOperatorCount(); i++) {

                    MCMCOperator op = schedule.getOperator(i);
                    double acceptanceProb = MCMCOperator.Utils.getAcceptanceProbability(op);
                    String message = "good";
                    if (acceptanceProb < op.getMinimumGoodAcceptanceLevel()) {
                        if (acceptanceProb < (op.getMinimumAcceptanceLevel()/10.0)) {
                            message = "very low";
                        } else if (acceptanceProb < op.getMinimumAcceptanceLevel()) {
                            message = "low";
                        } else message = "slightly low";

                    } else if (acceptanceProb > op.getMaximumGoodAcceptanceLevel()) {
                        double reallyHigh = 1.0 - ((1.0-op.getMaximumAcceptanceLevel())/10.0);
                        if (acceptanceProb > reallyHigh) {
                            message = "very high";
                        } else if (acceptanceProb > op.getMaximumAcceptanceLevel()) {
                            message = "high";
                        } else message = "slightly high";
                    }

                    String suggestion = op.getPerformanceSuggestion();

                    String pString = "        ";
                    if (op instanceof CoercableMCMCOperator && ((CoercableMCMCOperator)op).getMode() != CoercableMCMCOperator.COERCION_OFF) {
                        pString = formatter.formatToFieldWidth(formatter.formatDecimal(((CoercableMCMCOperator)op).getRawParameter(), 3), 8);
                    }

                    System.out.println(
                            formatter.formatToFieldWidth(op.getOperatorName(), 30) +

                                    pString +

                                    formatter.formatToFieldWidth(formatter.formatDecimal(acceptanceProb, 4), 11) +
                                    " " + message + "\t" + suggestion);
                }
                System.out.println();
            }

            // How should premature finish be flagged?
        }

    };

    /** @return the prior of this MCMC analysis. */
    public Prior getPrior() { return mc.getPrior(); }

    /** @return the likelihood function. */
    public Likelihood getLikelihood() { return mc.getLikelihood(); }

    /** @return the timer. */
    public dr.util.Timer getTimer() { return timer; }

    /** @return the length of this analysis.*/
    public final int getChainLength() { return options.getChainLength(); }

    // TRANSIENT PUBLIC METHODS *****************************************

    /** @return the current state of the MCMC analysis. */
    public final int getCurrentState() { return currentState; }

    /** @return the progress (0 to 1) of the MCMC analysis. */
    public final double getProgress() {
        return (double)currentState / (double)options.getChainLength();
    }

    /** @return true if this MCMC is currently adapting the operators. */
    public final boolean isAdapting() { return isAdapting; }

    /** Requests that the MCMC chain stop prematurely. */
    public void pleaseStop() {
        stopping = true;
        mc.pleaseStop();
    }

    /** sets the pause state of this MCMC chain */
    public boolean isStopped() { return mc.isStopped(); }

    //PRIVATE METHODS *****************************************
    private boolean isPreBurninNeeded() {

        if (options.useCoercion()) return true;

        for (int i = 0; i < schedule.getOperatorCount(); i++) {
            MCMCOperator op = schedule.getOperator(i);

            if (op instanceof CoercableMCMCOperator) {
                if (((CoercableMCMCOperator)op).getMode() == CoercableMCMCOperator.COERCION_ON) return true;
            }
        }
        return false;
    }

    public void setShowOperatorAnalysis(boolean soa) {
        showOperatorAnalysis = soa;
    }

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public String getParserName() { return MCMC; }

        /**
         * @return a tree object based on the XML element it was passed.
         */
        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            MCMC mcmc = new MCMC("mcmc1");
            MCMCOptions options = new MCMCOptions();
            OperatorSchedule opsched = (OperatorSchedule)xo.getChild(OperatorSchedule.class);
            Likelihood likelihood = (Likelihood)xo.getChild(Likelihood.class);
            ArrayList loggers = new ArrayList();

            options.setChainLength(xo.getIntegerAttribute(CHAIN_LENGTH));

            if (xo.hasAttribute(COERCION)) {
                options.setUseCoercion(xo.getBooleanAttribute(COERCION));
            }

            if (xo.hasAttribute(PRE_BURNIN)) {
                options.setPreBurnin(xo.getIntegerAttribute(PRE_BURNIN));
            }

            if (xo.hasAttribute(TEMPERATURE)) {
                options.setTemperature(xo.getIntegerAttribute(TEMPERATURE));
            }

            for (int i = 0; i < xo.getChildCount(); i++) {
                Object child = xo.getChild(i);
                if (child instanceof Logger) {
                    loggers.add((Logger)child);
                }
            }

            mcmc.setShowOperatorAnalysis(true);

            Logger[] loggerArray = new Logger[loggers.size()];
            loggers.toArray(loggerArray);

            java.util.logging.Logger.getLogger("dr.inference").info("Creating the MCMC chain:" +
                    "\n  chainLength=" + options.getChainLength() +
                    "\n  autoOptimize=" + options.useCoercion());

            mcmc.init(options, likelihood, Prior.UNIFORM_PRIOR, opsched, loggerArray);

            MarkovChain mc = mcmc.getMarkovChain();
            double initialScore = mc.getCurrentScore();

            if (initialScore == Double.NEGATIVE_INFINITY) {
                if (likelihood instanceof CompoundLikelihood) {
                    String message = ((CompoundLikelihood)likelihood).getDiagnosis();

                    throw new IllegalArgumentException("The initial posterior is zero: " + message);
                } else {
                    throw new IllegalArgumentException("The initial posterior is zero!");
                }
            }


            return mcmc;
        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public String getParserDescription() {
            return "This element returns an MCMC chain and runs the chain as a side effect.";
        }

        public Class getReturnType() { return MCMC.class; }

        public XMLSyntaxRule[] getSyntaxRules() { return rules; }

        private XMLSyntaxRule[] rules = new XMLSyntaxRule[] {
                AttributeRule.newIntegerRule(CHAIN_LENGTH),
                AttributeRule.newBooleanRule(COERCION, true),
                AttributeRule.newIntegerRule(PRE_BURNIN, true),
                AttributeRule.newDoubleRule(TEMPERATURE, true),
                new ElementRule(OperatorSchedule.class ),
                new ElementRule(Likelihood.class ),
                new ElementRule(Logger.class, 1, Integer.MAX_VALUE )
        };

    };

    public String getId() { return id; }

    public void setId(String id) { this.id = id; }

    // PRIVATE TRANSIENTS

    //private FileLogger operatorLogger = null;
    private boolean isAdapting = true, stopping = false;
    private boolean showOperatorAnalysis = true;
    private dr.util.Timer timer = new dr.util.Timer();
    private int currentState = 0;
    private int stepsPerReport = 1000;
    private NumberFormatter formatter = new NumberFormatter(8);

    /** this markov chain does most of the work. */
    private MarkovChain mc;

    /** the options of this MCMC analysis */
    private MCMCOptions options;

    private Logger[] loggers;
    private OperatorSchedule schedule;

    private String id = null;

    public static final String COERCION = "autoOptimize";
    public static final String PRE_BURNIN = "preBurnin";
    public static final String MCMC = "mcmc";
    public static final String CHAIN_LENGTH = "chainLength";
    public static final String WEIGHT = "weight";
    public static final String TEMPERATURE = "temperature";
}

