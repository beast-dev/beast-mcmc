/*
 * MCMCParser.java
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

package dr.inferencexml;

import dr.inference.loggers.Logger;
import dr.inference.markovchain.MarkovChain;
import dr.inference.mcmc.MCMC;
import dr.inference.mcmc.MCMCOptions;
import dr.inference.model.CompoundLikelihood;
import dr.inference.model.Likelihood;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.operators.OperatorSchedule;
import dr.xml.*;

import java.util.ArrayList;

public class MCMCParser extends AbstractXMLObjectParser {

    public String getParserName() {
        return MCMC;
    }

    /**
     * @return an mcmc object based on the XML element it was passed.
     */
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        MCMC mcmc = new MCMC(xo.getAttribute(NAME, "mcmc1"));

        long chainLength = xo.getLongIntegerAttribute(CHAIN_LENGTH);
        boolean useCoercion = xo.getAttribute(COERCION, true);
        long coercionDelay = chainLength / 100;
        if (xo.hasAttribute(PRE_BURNIN)) {
            coercionDelay = xo.getIntegerAttribute(PRE_BURNIN);
        }
        coercionDelay = xo.getAttribute(COERCION_DELAY, coercionDelay);
        double temperature = xo.getAttribute(TEMPERATURE, 1.0);

        long fullEvaluationCount = 1000;
        if (System.getProperty("mcmc.evaluation.count") != null) {
            fullEvaluationCount = Long.parseLong(System.getProperty("mcmc.evaluation.count"));
        }
        fullEvaluationCount = xo.getAttribute(FULL_EVALUATION, fullEvaluationCount);

        double evaluationTestThreshold = MarkovChain.EVALUATION_TEST_THRESHOLD;
        if (System.getProperty("mcmc.evaluation.threshold") != null) {
            evaluationTestThreshold = Double.parseDouble(System.getProperty("mcmc.evaluation.threshold"));
        }
        evaluationTestThreshold = xo.getAttribute(EVALUATION_THRESHOLD, evaluationTestThreshold);

        int minOperatorCountForFullEvaluation = xo.getAttribute(MIN_OPS_EVALUATIONS, 1);

        MCMCOptions options = new MCMCOptions(chainLength,
                fullEvaluationCount,
                minOperatorCountForFullEvaluation,
                evaluationTestThreshold,
                useCoercion,
                coercionDelay,
                temperature);

        OperatorSchedule opsched = (OperatorSchedule) xo.getChild(OperatorSchedule.class);
        Likelihood likelihood = (Likelihood) xo.getChild(Likelihood.class);

        likelihood.setUsed();

        if (Boolean.valueOf(System.getProperty("show_warnings", "false"))) {

            // check that all models, parameters and likelihoods are being used
            for (Likelihood l : Likelihood.FULL_LIKELIHOOD_SET) {
                if (!l.isUsed()) {
                    java.util.logging.Logger.getLogger("dr.inference").warning("Likelihood, " + l.getId() +
                            ", of class " + l.getClass().getName() + " is not being handled by the MCMC.");
                }
            }
            for (Model m : Model.FULL_MODEL_SET) {
                if (!m.isUsed()) {
                    java.util.logging.Logger.getLogger("dr.inference").warning("Model, " + m.getId() +
                            ", of class " + m.getClass().getName() + " is not being handled by the MCMC.");
                }
            }
            for (Parameter p : Parameter.FULL_PARAMETER_SET) {
                if (!p.isUsed()) {
                    java.util.logging.Logger.getLogger("dr.inference").warning("Parameter, " + p.getId() +
                            ", of class " + p.getClass().getName() + " is not being handled by the MCMC.");
                }
            }
        }

        ArrayList<Logger> loggers = new ArrayList<Logger>();

        for (int i = 0; i < xo.getChildCount(); i++) {
            Object child = xo.getChild(i);
            if (child instanceof Logger) {
                loggers.add((Logger) child);
            }
        }

        mcmc.setShowOperatorAnalysis(true);
        if (xo.hasAttribute(OPERATOR_ANALYSIS)) {
            mcmc.setOperatorAnalysisFile(XMLParser.getLogFile(xo, OPERATOR_ANALYSIS));
        }


        Logger[] loggerArray = new Logger[loggers.size()];
        loggers.toArray(loggerArray);


        java.util.logging.Logger.getLogger("dr.inference").info("\nCreating the MCMC chain:" +
                "\n  chainLength=" + options.getChainLength() +
                "\n  autoOptimize=" + options.useCoercion() +
                (options.useCoercion() ? "\n  autoOptimize delayed for " + options.getCoercionDelay() + " steps" : "") +
                (options.getFullEvaluationCount() == 0 ? "\n  full evaluation test off" : "")
        );

        mcmc.init(options, likelihood, opsched, loggerArray);


        MarkovChain mc = mcmc.getMarkovChain();
        double initialScore = mc.getCurrentScore();

        if (initialScore == Double.NEGATIVE_INFINITY) {
            String message = "The initial posterior is zero";
            if (likelihood instanceof CompoundLikelihood) {
                message += ": " + ((CompoundLikelihood) likelihood).getDiagnosis(2);
            } else {
                message += "!";
            }
            throw new IllegalArgumentException(message);
        }

        if (!xo.getAttribute(SPAWN, true))
            mcmc.setSpawnable(false);

        return mcmc;
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "This element returns an MCMC chain and runs the chain as a side effect.";
    }

    public Class getReturnType() {
        return MCMC.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            AttributeRule.newLongIntegerRule(CHAIN_LENGTH),
            AttributeRule.newBooleanRule(COERCION, true),
            AttributeRule.newIntegerRule(COERCION_DELAY, true),
            AttributeRule.newIntegerRule(PRE_BURNIN, true),
            AttributeRule.newDoubleRule(TEMPERATURE, true),
            AttributeRule.newIntegerRule(FULL_EVALUATION, true),
            AttributeRule.newIntegerRule(MIN_OPS_EVALUATIONS, true),
            AttributeRule.newDoubleRule(EVALUATION_THRESHOLD, true),
            AttributeRule.newBooleanRule(SPAWN, true),
            AttributeRule.newStringRule(NAME, true),
            AttributeRule.newStringRule(OPERATOR_ANALYSIS, true),
            new ElementRule(OperatorSchedule.class),
            new ElementRule(Likelihood.class),
            new ElementRule(Logger.class, 1, Integer.MAX_VALUE),
    };

    public static final String COERCION = "autoOptimize";
    public static final String NAME = "name";
    public static final String PRE_BURNIN = "preBurnin";
    public static final String COERCION_DELAY = "autoOptimizeDelay";
    public static final String MCMC = "mcmc";
    public static final String CHAIN_LENGTH = "chainLength";
    public static final String FULL_EVALUATION = "fullEvaluation";
    public static final String EVALUATION_THRESHOLD  = "evaluationThreshold";
    public static final String MIN_OPS_EVALUATIONS = "minOpsFullEvaluations";
    public static final String WEIGHT = "weight";
    public static final String TEMPERATURE = "temperature";
    public static final String SPAWN = "spawn";
    public static final String OPERATOR_ANALYSIS = "operatorAnalysis";


}
