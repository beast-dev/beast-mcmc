/*
 * MarkovChain.java
 *
 * Copyright (c) 2002-2012 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.inference.markovchain;

import dr.inference.model.CompoundLikelihood;
import dr.inference.model.Likelihood;
import dr.inference.model.Model;
import dr.inference.operators.*;
import dr.inference.prior.Prior;

import java.util.ArrayList;
import java.util.logging.Logger;

/**
 * A concrete markov chain. This is final as the only things that should need
 * overriding are in the delegates (prior, likelihood, schedule and acceptor).
 * The design of this class is to be fairly immutable as far as settings goes.
 *
 * @author Alexei Drummond
 * @author Andrew Rambaut
 * @version $Id: MarkovChain.java,v 1.10 2006/06/21 13:34:42 rambaut Exp $
 */
public final class MarkovChain {

    private final static boolean DEBUG = false;
    private final static boolean PROFILE = true;

    private final OperatorSchedule schedule;
    private final Acceptor acceptor;
    private final Prior prior;
    private final Likelihood likelihood;

    private boolean pleaseStop = false;
    private boolean isStopped = false;
    private double bestScore, currentScore, initialScore;
    private long currentLength;

    private boolean useCoercion = true;

    private final int fullEvaluationCount;
    private final int minOperatorCountForFullEvaluation;

    private static final double EVALUATION_TEST_THRESHOLD = 1e-6;

    public MarkovChain(Prior prior, Likelihood likelihood,
                       OperatorSchedule schedule, Acceptor acceptor,
                       int fullEvaluationCount, int minOperatorCountForFullEvaluation, boolean useCoercion) {
        currentLength = 0;
        this.prior = prior;
        this.likelihood = likelihood;
        this.schedule = schedule;
        this.acceptor = acceptor;
        this.useCoercion = useCoercion;

        this.fullEvaluationCount = fullEvaluationCount;
        this.minOperatorCountForFullEvaluation = minOperatorCountForFullEvaluation;

        currentScore = evaluate(likelihood, prior);
    }

    /**
     * Resets the markov chain
     */
    public void reset() {
        currentLength = 0;

        // reset operator acceptance levels
        for (int i = 0; i < schedule.getOperatorCount(); i++) {
            schedule.getOperator(i).reset();
        }
    }

    /**
     * Run the chain for a given number of states.
     *
     * @param length number of states to run the chain.
     *               <p/>
     *               param onTheFlyOperatorWeights
     */
    public long runChain(long length, boolean disableCoerce /*,int onTheFlyOperatorWeights*/) {

        likelihood.makeDirty();
        currentScore = evaluate(likelihood, prior);

        long currentState = currentLength;

        final Model currentModel = likelihood.getModel();

        if (currentState == 0) {
            initialScore = currentScore;
            bestScore = currentScore;
            fireBestModel(currentState, currentModel);
        }

        if (currentScore == Double.NEGATIVE_INFINITY) {

            // identify which component of the score is zero...
            if (prior != null) {
                double logPrior = prior.getLogPrior(likelihood.getModel());

                if (logPrior == Double.NEGATIVE_INFINITY) {
                    throw new IllegalArgumentException(
                            "The initial model is invalid because one of the priors has zero probability.");
                }
            }

            String message = "The initial likelihood is zero";
            if (likelihood instanceof CompoundLikelihood) {
                message += ": " + ((CompoundLikelihood) likelihood).getDiagnosis();
            } else {
                message += ".";
            }
            throw new IllegalArgumentException(message);
        } else if (currentScore == Double.POSITIVE_INFINITY || Double.isNaN(currentScore)) {
            String message = "A likelihood returned with a numerical error";
            if (likelihood instanceof CompoundLikelihood) {
                message += ": " + ((CompoundLikelihood) likelihood).getDiagnosis();
            } else {
                message += ".";
            }
            throw new IllegalArgumentException(message);
        }

        pleaseStop = false;
        isStopped = false;

        String diagnostic = "";

        //int otfcounter = onTheFlyOperatorWeights > 0 ? onTheFlyOperatorWeights : 0;

        double[] logr = {0.0};

        boolean usingFullEvaluation = true;
        // set ops count in mcmc element instead
        if (fullEvaluationCount == 0) // Temporary solution until full code review
            usingFullEvaluation = false;
        boolean fullEvaluationError = false;

        while (!pleaseStop && (currentState < (currentLength + length))) {

            // periodically log states
            fireCurrentModel(currentState, currentModel);

            if (pleaseStop) {
                isStopped = true;
                break;
            }

            // Get the operator
            final int op = schedule.getNextOperatorIndex();
            final MCMCOperator mcmcOperator = schedule.getOperator(op);

            double oldScore = currentScore;

            // assert Profiler.startProfile("Store");

            // The current model is stored here in case the proposal fails
            if (currentModel != null) {
                currentModel.storeModelState();
            }

            // assert Profiler.stopProfile("Store");

            boolean operatorSucceeded = true;
            double hastingsRatio = 1.0;
            boolean accept = false;

            logr[0] = -Double.MAX_VALUE;

            try {
                // The new model is proposed
                // assert Profiler.startProfile("Operate");

                if (DEBUG) {
                    System.out.println("\n&& Operator: " + mcmcOperator.getOperatorName());
                }

                if (mcmcOperator instanceof GeneralOperator) {
                    hastingsRatio = ((GeneralOperator) mcmcOperator).operate(prior, likelihood);
                } else {
                    hastingsRatio = mcmcOperator.operate();
                }

                // assert Profiler.stopProfile("Operate");
            } catch (OperatorFailedException e) {
                operatorSucceeded = false;
            }

            double score = 0.0;
            double deviation = 0.0;

            //    System.err.print("" + currentState + ": ");
            if (operatorSucceeded) {

                // The new model is proposed
                // assert Profiler.startProfile("Evaluate");

                if (DEBUG) {
                    System.out.println("** Evaluate");
                }

                long elapsedTime = 0;
                if (PROFILE) {
                    elapsedTime = System.currentTimeMillis();
                }

                // The new model is evaluated
                score = evaluate(likelihood, prior);

                if (PROFILE) {
                    mcmcOperator.addEvaluationTime(System.currentTimeMillis() - elapsedTime);
                }

                // assert Profiler.stopProfile("Evaluate");

                if (score == Double.POSITIVE_INFINITY || Double.isNaN(score)) {
                    if (likelihood instanceof CompoundLikelihood) {
                        Logger.getLogger("error").severe("A likelihood returned with a numerical error:\n" +
                                ((CompoundLikelihood)likelihood).getDiagnosis());
                    } else {
                        Logger.getLogger("error").severe("A likelihood returned with a numerical error.");
                    }

                    // If the user has chosen to ignore this error then we transform it
                    // to a negative infinity so the state is rejected.
                    score = Double.NEGATIVE_INFINITY;
                }

                if (usingFullEvaluation) {

                    // This is a test that the state is correctly restored. The
                    // restored state is fully evaluated and the likelihood compared with
                    // that before the operation was made.
                    likelihood.makeDirty();
                    final double testScore = evaluate(likelihood, prior);

                    if (Math.abs(testScore - score) > EVALUATION_TEST_THRESHOLD) {
                        Logger.getLogger("error").severe(
                                "State was not correctly calculated after an operator move.\n"
                                        + "Likelihood evaluation: " + score
                                        + "\nFull Likelihood evaluation: " + testScore
                                        + "\n" + "Operator: " + mcmcOperator
                                        + " " + mcmcOperator.getOperatorName());
                        fullEvaluationError = true;
                    }
                }

                if (score > bestScore) {
                    bestScore = score;
                    fireBestModel(currentState, currentModel);
                }

                accept = mcmcOperator instanceof GibbsOperator || acceptor.accept(oldScore, score, hastingsRatio, logr);

                deviation = score - oldScore;
            }

            // The new model is accepted or rejected
            if (accept) {
                if (DEBUG) {
                    System.out.println("** Move accepted: new score = " + score
                            + ", old score = " + oldScore);
                }

                mcmcOperator.accept(deviation);
                currentModel.acceptModelState();
                currentScore = score;

//                if( otfcounter > 0 ) {
//                    --otfcounter;
//                    if( otfcounter == 0 ) {
//                        adjustOpWeights(currentState);
//                        otfcounter = onTheFlyOperatorWeights;
//                    }
//                }

                if (usingFullEvaluation) {
                    oldScore = score; // for the usingFullEvaluation test
                    diagnostic = likelihood instanceof CompoundLikelihood ?
                            ((CompoundLikelihood) likelihood).getDiagnosis() : "";
                }
            } else {
                if (DEBUG) {
                    System.out.println("** Move rejected: new score = " + score
                            + ", old score = " + oldScore);
                }

                mcmcOperator.reject();

                // assert Profiler.startProfile("Restore");

                currentModel.restoreModelState();
            }
            // assert Profiler.stopProfile("Restore");

            if (usingFullEvaluation) {
                // This is a test that the state is correctly restored. The
                // restored state is fully evaluated and the likelihood compared with
                // that before the operation was made.

                likelihood.makeDirty();
                final double testScore = evaluate(likelihood, prior);

                final String d2 = likelihood instanceof CompoundLikelihood ?
                        ((CompoundLikelihood) likelihood).getDiagnosis() : "";

                if (Math.abs(testScore - oldScore) > EVALUATION_TEST_THRESHOLD) {


                    final Logger logger = Logger.getLogger("error");
                    logger.severe("State was not correctly restored after reject step.\n"
                            + "Likelihood before: " + oldScore
                            + " Likelihood after: " + testScore
                            + "\n" + "Operator: " + mcmcOperator
                            + " " + mcmcOperator.getOperatorName()
                            + (diagnostic.length() > 0 ? "\n\nDetails\nBefore: " + diagnostic + "\nAfter: " + d2 : "")
                    );
                    fullEvaluationError = true;
                }
            }

            if (!disableCoerce && mcmcOperator instanceof CoercableMCMCOperator) {
                coerceAcceptanceProbability((CoercableMCMCOperator) mcmcOperator, logr[0]);
            }

            if (usingFullEvaluation) {
                if (schedule.getMinimumAcceptAndRejectCount() >= minOperatorCountForFullEvaluation &&
                        currentState >= fullEvaluationCount) {
                    // full evaluation is only switched off when each operator has done a
                    // minimum number of operations (currently 1) and fullEvalationCount
                    // operations in total.

                    usingFullEvaluation = false;
                    if (fullEvaluationError) {
                        // If there has been an error then stop with an error
                        throw new RuntimeException(
                                "One or more evaluation errors occurred during the test phase of this\n" +
                                        "run. These errors imply critical errors which may produce incorrect\n" +
                                        "results.");
                    }
                }
            }

            fireEndCurrentIteration(currentState);

            currentState += 1;
        }

        currentLength = currentState;

        return currentLength;
    }

    public void terminateChain() {
        fireFinished(currentLength);

        // Profiler.report();
    }
//
//    private void adjustOpWeights(int currentState) {
//        final int count = schedule.getOperatorCount();
//        double[] s = new double[count];
//        final double factor = 100;
//        final double limitSpan = 1000;
//        System.err.println("start cycle " + currentState);
//
//        double sHas = 0.0/* , sNot = 0.0 */, nHas = 0.0;
//        for(int no = 0; no < count; ++no) {
//            final MCMCOperator op = schedule.getOperator(no);
//            final double v = op.getSpan(true);
//
//            if( v == 0 ) {
//                // sNot += op.getWeight();
//                s[no] = 0;
//            } else {
//                sHas += op.getWeight();
//                s[no] = Math.max(factor * Math.min(v, limitSpan), 1);
//                nHas += s[no];
//            }
//        }
//
//        // for(int no = 0; no < count; ++no) {
//        // final MCMCOperator op = schedule.getOperator(no);
//        // final double v = op.getSpan(false);
//        // if( v == 0 ) {
//        // System.err.println(op.getOperatorName() + " blocks");
//        // return;
//        // }
//        // }
//
//        // keep sum of changed parts unchanged
//        final double scaleHas = sHas / nHas;
//
//        for(int no = 0; no < count; ++no) {
//            final MCMCOperator op = schedule.getOperator(no);
//            if( s[no] > 0 ) {
//                final double val = s[no] * scaleHas;
//                op.setWeight(val);
//                System.err.println("set " + op.getOperatorName() + " " + val);
//            } else {
//                System.err.println("** " + op.getOperatorName() + " = "
//                        + op.getWeight());
//            }
//        }
//        schedule.operatorsHasBeenUpdated();
//    }


    public Prior getPrior() {
        return prior;
    }

    public Likelihood getLikelihood() {
        return likelihood;
    }

    public Model getModel() {
        return likelihood.getModel();
    }

    public OperatorSchedule getSchedule() {
        return schedule;
    }

    public Acceptor getAcceptor() {
        return acceptor;
    }

    public double getInitialScore() {
        return initialScore;
    }

    public double getBestScore() {
        return bestScore;
    }

    public long getCurrentLength() {
        return currentLength;
    }

    public void setCurrentLength(long currentLength) {
        this.currentLength = currentLength;
    }

    public double getCurrentScore() {
        return currentScore;
    }

    public void pleaseStop() {
        pleaseStop = true;
    }

    public boolean isStopped() {
        return isStopped;
    }

    protected double evaluate(Likelihood likelihood, Prior prior) {

        double logPosterior = 0.0;

        if (prior != null) {
            final double logPrior = prior.getLogPrior(likelihood.getModel());

            if (logPrior == Double.NEGATIVE_INFINITY) {
                return Double.NEGATIVE_INFINITY;
            }

            logPosterior += logPrior;
        }

        final double logLikelihood = likelihood.getLogLikelihood();

        if (Double.isNaN(logLikelihood)) {
            return Double.NEGATIVE_INFINITY;
        }
        // System.err.println("** " + logPosterior + " + " + logLikelihood +
        // " = " + (logPosterior + logLikelihood));
        logPosterior += logLikelihood;

        return logPosterior;
    }

    /**
     * Updates the proposal parameter, based on the target acceptance
     * probability This method relies on the proposal parameter being a
     * decreasing function of acceptance probability.
     *
     * @param op   The operator
     * @param logr
     */
    private void coerceAcceptanceProbability(CoercableMCMCOperator op, double logr) {

        if (isCoercable(op)) {
            final double p = op.getCoercableParameter();

            final double i = schedule.getOptimizationTransform(MCMCOperator.Utils.getOperationCount(op));

            final double target = op.getTargetAcceptanceProbability();

            final double newp = p + ((1.0 / (i + 1.0)) * (Math.exp(logr) - target));

            if (newp > -Double.MAX_VALUE && newp < Double.MAX_VALUE) {
                op.setCoercableParameter(newp);
            }
        }
    }

    private boolean isCoercable(CoercableMCMCOperator op) {

        return op.getMode() == CoercionMode.COERCION_ON
                || (op.getMode() != CoercionMode.COERCION_OFF && useCoercion);
    }

    public void addMarkovChainListener(MarkovChainListener listener) {
        listeners.add(listener);
    }

    public void removeMarkovChainListener(MarkovChainListener listener) {
        listeners.remove(listener);
    }

    public void addMarkovChainDelegate(MarkovChainDelegate delegate) {
        delegates.add(delegate);
    }

    public void removeMarkovChainDelegate(MarkovChainDelegate delegate) {
        delegates.remove(delegate);
    }


    private void fireBestModel(long state, Model bestModel) {

        for (MarkovChainListener listener : listeners) {
            listener.bestState(state, bestModel);
        }
    }

    private void fireCurrentModel(long state, Model currentModel) {
        for (MarkovChainListener listener : listeners) {
            listener.currentState(state, currentModel);
        }

        for (MarkovChainDelegate delegate : delegates) {
            delegate.currentState(state);
        }
    }

    private void fireFinished(long chainLength) {

        for (MarkovChainListener listener : listeners) {
            listener.finished(chainLength);
        }

        for (MarkovChainDelegate delegate : delegates) {
            delegate.finished(chainLength);
        }
    }

    private void fireEndCurrentIteration(long state) {
        for (MarkovChainDelegate delegate : delegates) {
            delegate.currentStateEnd(state);
        }
    }

    private final ArrayList<MarkovChainListener> listeners = new ArrayList<MarkovChainListener>();
    private final ArrayList<MarkovChainDelegate> delegates = new ArrayList<MarkovChainDelegate>();
}
