/*
 * BlueBeastMCMCParser.java
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

/*
 * Interaction with BEAST XML
 * XML options to turn our program on
 * Warning message when turned on: Please cite our paper if you use this. Advise to check chain manually
 * Weights placed on the proposal kernels will be treated as starting weights
 * Optional XML options (hidden) -
 *      constant/dynamic interval change
 *      autooptimise weights
 *      maximum chain length (if not infinity)
 *      convergence assessment statistics (default is to use all)
 *      load tracer or not when convergence is reached
 *
 *
 * Mostly copied from MCMCParser.java
 *
 * @author Wai Lok Sibon Li

*/

package dr.inferencexml;
//package beast.parser;

import bb.main.BlueBeast;
import bb.mcmc.analysis.*;
import beast.core.BlueBeastMCMC;
import beast.core.BlueBeastMarkovChain;
import beast.inference.loggers.BlueBeastLogger;
import dr.inference.loggers.Logger;
import dr.inference.loggers.MCLogger;
import dr.inference.mcmc.MCMC;
import dr.inference.mcmc.MCMCOptions;
import dr.inference.model.CompoundLikelihood;
import dr.inference.model.Likelihood;
import dr.inference.operators.OperatorSchedule;
import dr.xml.*;

import java.util.ArrayList;


/**
 * This class is actually not used anywhere within BlueBEAST, instead, a copy of this class can be found in the
 * BEAST source code in dr.inferencexml
 */
public class BlueBeastMCMCParser extends AbstractXMLObjectParser {

    public String getParserName() {
        return MCMC;
    }

    /**
     * @return an mcmc object based on the XML element it was passed.
     */
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        BlueBeastMCMC mcmc = new BlueBeastMCMC(xo.getAttribute(NAME, "mcmc1"));
        MCMCOptions options = new MCMCOptions();
        OperatorSchedule opsched = (OperatorSchedule) xo.getChild(OperatorSchedule.class);
        Likelihood likelihood = (Likelihood) xo.getChild(Likelihood.class);
//        BlueBeast bb = (BlueBeast) xo.getChild(BlueBeast.class);
        BlueBeastLogger bbl = (BlueBeastLogger) xo.getChild(BlueBeastLogger.class);

        ArrayList<Logger> loggers = new ArrayList<Logger>();

        likelihood.setUsed();

        //options.setChainLength(xo.getIntegerAttribute(CHAIN_LENGTH));
        int maxChainLength = xo.getAttribute(MAX_CHAIN_LENGTH, Integer.MAX_VALUE);
        //options.setChainLength(maxChainLength);

        options.setUseCoercion(xo.getAttribute(COERCION, true));
        if (xo.hasAttribute(PRE_BURNIN)) {
            options.setCoercionDelay(xo.getAttribute(PRE_BURNIN, maxChainLength / 100));
        }
        options.setCoercionDelay(xo.getAttribute(COERCION_DELAY, maxChainLength / 100));
        options.setTemperature(xo.getAttribute(TEMPERATURE, 1.0));
        options.setFullEvaluationCount(xo.getAttribute(FULL_EVALUATION, 2000));
        options.setMinOperatorCountForFullEvaluation(xo.getAttribute(MIN_OPS_EVALUATIONS, 1));

//        if (xo.hasAttribute(CHECK_INTERVAL)) {
//            mcmc.setCheckInterval(xo.getIntegerAttribute(CHECK_INTERVAL));
//        }
//        else {
//            mcmc.setCheckInterval(1000);
//        }


        long initialCheckInterval = 1000;
        if (xo.hasAttribute(INITIAL_CHECK_INTERVAL)) {
            initialCheckInterval = xo.getIntegerAttribute(INITIAL_CHECK_INTERVAL);
            System.out.println("Starting log interval: " + initialCheckInterval);
        }
        mcmc.setCheckInterval(initialCheckInterval);

//        ArrayList<String> varNames = new ArrayList<String>(xo.getChildCount());
        for (int i = 0; i < xo.getChildCount(); i++) {
            Object child = xo.getChild(i);
            if (child instanceof Logger) {
                loggers.add((Logger) child);
                // TODO Get file names
                //Logger temp = (Logger) child;
                //if(temp.)
                //loggers.add(temp);
                //if(child instanceof MCLogger) {
                //    MCLogger mcLogger = (MCLogger) child;
                //    mcLogger.getTitle()
                //}

//                    for(int j = 0; j < mcLogger.getColumnCount(); j++) {
//                        varNames.add(mcLogger.getColumnLabel(j));
////                        System.out.println("test 4 " + mcLogger.getColumnLabel(j));
//                    }
////                    if(mcLogger.getFormatters().get(0) instanceof TabDelimitedFormatter) {
////                        System.out.println("check check");
////                    }
//
//                    //varNames.add(((MCLogger) child).getTitle()); // Hope this works
//                    //((MCLogger) child).getColumnLabel(0)
//                    //varNames.add(((MCLogger) child).getColumn(0).getLabel()); // Hope this works
//                }

            }
        }

        mcmc.setShowOperatorAnalysis(true);
        if (xo.hasAttribute(OPERATOR_ANALYSIS)) {
            mcmc.setOperatorAnalysisFileName(xo.getStringAttribute(OPERATOR_ANALYSIS));
        }

        Logger[] loggerArray = new Logger[loggers.size()];
        loggers.toArray(loggerArray);

        java.util.logging.Logger.getLogger("dr.inference").info("Creating the BlueBeastMCMC chain:" +
                "\n  maxChainLength=" + maxChainLength +
                "\n  autoOptimize=" + options.useCoercion() +
                (options.useCoercion() ? "\n  autoOptimize delayed for " + options.getCoercionDelay() + " steps" : "") +
                (options.fullEvaluationCount() == 0 ? "\n  full evaluation test off" : "")
                );


        ArrayList<ConvergeStat> convergenceStatsToUse = new ArrayList<ConvergeStat>();
        String convergenceStatsToUseParameters = xo.getAttribute(CONVERGENCE_STATS_TO_USE, "all");
        if(convergenceStatsToUseParameters.equals("all")) {
            convergenceStatsToUse.add(ESSConvergeStat.INSTANCE);
            convergenceStatsToUse.add(GelmanConvergeStat.INSTANCE);
            convergenceStatsToUse.add(GewekeConvergeStat.INSTANCE);
            convergenceStatsToUse.add(ZTempNovelConvergenceStatistic.INSTANCE);
        }
        else if(convergenceStatsToUseParameters.equals("ESS")) {
            convergenceStatsToUse.add(ESSConvergeStat.INSTANCE);
        }
        if(convergenceStatsToUseParameters.equals("interIntraChainVariance")) {
            convergenceStatsToUse.add(GewekeConvergeStat.INSTANCE);
            convergenceStatsToUse.add(GelmanConvergeStat.INSTANCE);
        }

        //String[] variableNames = varNames.toArray(new String[varNames.size()]);
//        String[] variableNames = null;//bbl.getvariableNames();

        int essLowerLimitBoundary = xo.getAttribute(ESS_LOWER_LIMIT_BOUNDARY, 200);
        double burninPercentage = xo.getAttribute(BURNIN_PERCENTAGE, 0.1);
        boolean dynamicCheckingInterval = xo.getAttribute(DYNAMIC_CHECKING_INTERVAL, true);
        boolean autoOptimiseWeights = xo.getAttribute(AUTO_OPTIMISE_WEIGHTS, true);
        boolean optimiseChainLength = xo.getAttribute(OPTIMISE_CHAIN_LENGTH, true);
        boolean loadTracer = xo.getAttribute(LOAD_TRACER, true);

        //options.setChainLength(xo.getIntegerAttribute(MAX_CHAIN_LENGTH));
        BlueBeast bb = new BlueBeast(opsched, options, convergenceStatsToUse, bbl, essLowerLimitBoundary,
                burninPercentage, dynamicCheckingInterval, autoOptimiseWeights, optimiseChainLength, maxChainLength,
                initialCheckInterval, loadTracer);

        mcmc.init(options, likelihood, opsched, loggerArray, bb);

        BlueBeastMarkovChain mc = (BlueBeastMarkovChain) mcmc.getMarkovChain();
        //MarkovChain mc = mcmc.getMarkovChain();
        if(mc == null) {
            System.out.println("shot bro");
        }
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
            //AttributeRule.newIntegerRule(CHAIN_LENGTH),
            AttributeRule.newBooleanRule(COERCION, true),
            AttributeRule.newIntegerRule(COERCION_DELAY, true),
            AttributeRule.newIntegerRule(PRE_BURNIN, true),
            AttributeRule.newDoubleRule(TEMPERATURE, true),
            AttributeRule.newIntegerRule(FULL_EVALUATION, true),
            AttributeRule.newIntegerRule(MIN_OPS_EVALUATIONS, true),
            AttributeRule.newBooleanRule(SPAWN, true),
            AttributeRule.newStringRule(NAME, true),
            AttributeRule.newStringRule(OPERATOR_ANALYSIS, true),



//            AttributeRule.newStringRule(CHECK_INTERVAL, true),
            AttributeRule.newLongIntegerRule(INITIAL_CHECK_INTERVAL, true),
            AttributeRule.newStringRule(CONVERGENCE_STATS_TO_USE, true),
            AttributeRule.newIntegerRule(ESS_LOWER_LIMIT_BOUNDARY, true),
            AttributeRule.newDoubleRule(BURNIN_PERCENTAGE, true),
            AttributeRule.newBooleanRule(DYNAMIC_CHECKING_INTERVAL, true),
            AttributeRule.newBooleanRule(AUTO_OPTIMISE_WEIGHTS, true),
            AttributeRule.newBooleanRule(OPTIMISE_CHAIN_LENGTH, true),
            AttributeRule.newBooleanRule(LOAD_TRACER, true),
            AttributeRule.newLongIntegerRule(MAX_CHAIN_LENGTH, true),

            new ElementRule(OperatorSchedule.class),
            new ElementRule(Likelihood.class),
            new ElementRule(BlueBeastLogger.class),
            new ElementRule(Logger.class, 1, Integer.MAX_VALUE),

//            new ElementRule(BlueBeast.class)
    };

    public static final String COERCION = "autoOptimize";
    public static final String NAME = "name";
    public static final String PRE_BURNIN = "preBurnin";
    public static final String COERCION_DELAY = "autoOptimizeDelay";
    public static final String MCMC = "blueBeastmcmc";
    //public static final String CHAIN_LENGTH = "chainLength";
    public static final String FULL_EVALUATION = "fullEvaluation";
    public static final String MIN_OPS_EVALUATIONS = "minOpsFullEvaluations";
    public static final String WEIGHT = "weight";
    public static final String TEMPERATURE = "temperature";
    public static final String SPAWN = "spawn";
    public static final String OPERATOR_ANALYSIS = "operatorAnalysis";
    //public static final String CHECK_INTERVAL = "checkInterval";
    public static final String CONVERGENCE_STATS_TO_USE = "convergenceStatsToUse";
    public static final String ESS_LOWER_LIMIT_BOUNDARY = "essLowerLimitBoundary";
    public static final String BURNIN_PERCENTAGE = "burninPercentage";
    public static final String DYNAMIC_CHECKING_INTERVAL = "dynamicCheckingInterval";
    public static final String AUTO_OPTIMISE_WEIGHTS = "autoOptimiseWeights";
    public static final String OPTIMISE_CHAIN_LENGTH = "optimiseChainLength";
    public static final String MAX_CHAIN_LENGTH = "maxChainLength";
    public static final String INITIAL_CHECK_INTERVAL = "initialCheckInterval";
    public static final String LOAD_TRACER = "loadTracer";







}
