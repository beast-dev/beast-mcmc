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
import dr.inference.loggers.MCLogger;
import dr.inference.markovchain.MarkovChain;
import dr.inference.markovchain.MarkovChainListener;
import dr.inference.model.Model;
import dr.inference.model.PathLikelihood;
import dr.inference.operators.OperatorSchedule;
import dr.inference.prior.Prior;
import dr.util.Identifiable;
import dr.xml.*;

import java.util.ArrayList;

/**
 * An MCMC analysis that estimates parameters of a probabilistic model.
 *
 * @author Alexei Drummond
 * @author Andrew Rambaut
 * @version $Id: MCMC.java,v 1.41 2005/07/11 14:06:25 rambaut Exp $
 */
public class MarginalLikelihoodEstimator implements Runnable, Identifiable {

    public MarginalLikelihoodEstimator(String id, int chainLength, int pathSteps, PathLikelihood pathLikelihood,
                                       OperatorSchedule schedule,
                                       MCLogger logger) {

        this.id = id;
        this.chainLength = chainLength;
        this.pathSteps = pathSteps;


        MCMCCriterion criterion = new MCMCCriterion();

        //pathDelta = 1.0 / pathSteps;
        pathParameter = 1.0;

        this.pathLikelihood = pathLikelihood;
        pathLikelihood.setPathParameter(pathParameter);

        mc = new MarkovChain(Prior.UNIFORM_PRIOR, pathLikelihood, schedule, criterion, 0, false);

        this.logger = logger;
    }

    public void run() {

        int logCount = chainLength / logger.getLogEvery();

        logger.startLogging();

        mc.addMarkovChainListener(chainListener);

        for (int step = 0; step < pathSteps; step++) {
            double p = 1.0 - (((double)step) / pathSteps);
            int cl = (int)(20.0 * chainLength / ((19.0 * p) + 1));
            logger.setLogEvery(cl / logCount);
            burnin = (int)(0.1 * cl);  

            pathLikelihood.setPathParameter(pathParameter);
            mc.setCurrentLength(0);
            mc.chain(cl + burnin, false, 0, false);
            pathParameter /= 2;

        }

        int cl = 20 * chainLength;
        logger.setLogEvery(cl / logCount);
        burnin = (int)(0.1 * cl);
        pathLikelihood.setPathParameter(0.0);
        mc.setCurrentLength(0);
        mc.chain(cl + burnin, false, 0, false);

        mc.removeMarkovChainListener(chainListener);
    }

    private MarkovChainListener chainListener = new MarkovChainListener() {

        // MarkovChainListener interface *******************************************
        // for receiving messages from subordinate MarkovChain

        /**
         * Called to update the current model keepEvery states.
         */
        public void currentState(int state, Model currentModel) {

            currentState = state;

            if (currentState >= burnin) {
            logger.log(state);
            }
        }

        /**
         * Called when a new new best posterior state is found.
         */
        public void bestState(int state, Model bestModel) {
            currentState = state;
        }

        /**
         * cleans up when the chain finishes (possibly early).
         */
        public void finished(int chainLength) {
            currentState = chainLength;

//            logger.log(currentState);
            logger.stopLogging();
        }
    };

    // TRANSIENT PUBLIC METHODS *****************************************

    /**
     * @return the current state of the MCMC analysis.
     */
    public boolean getSpawnable() {
        return spawnable;
    }

    private boolean spawnable = true;

    public void setSpawnable(boolean spawnable) {
        this.spawnable = spawnable;
    }


    //PRIVATE METHODS *****************************************
    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return MARGINAL_LIKELIHOOD_ESTIMATOR;
        }

        /**
         * @return a tree object based on the XML element it was passed.
         */
        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            MCMC mcmc = (MCMC) xo.getChild(MCMC.class);

            PathLikelihood pathLikelihood = (PathLikelihood) xo.getChild(PathLikelihood.class);
            MCLogger logger = (MCLogger) xo.getChild(MCLogger.class);

            int chainLength = xo.getIntegerAttribute(CHAIN_LENGTH);
            int pathSteps = xo.getIntegerAttribute(PATH_STEPS);

            for (int i = 0; i < xo.getChildCount(); i++) {
                Object child = xo.getChild(i);
                if (child instanceof Logger) {
                }
            }

            java.util.logging.Logger.getLogger("dr.inference").info("Creating the Marginal Likelihood Estimator chain:" +
                    "\n  chainLength=" + chainLength +
                    "\n  pathSteps=" + pathSteps);

            MarginalLikelihoodEstimator mle = new MarginalLikelihoodEstimator(MARGINAL_LIKELIHOOD_ESTIMATOR, chainLength,
                    pathSteps, pathLikelihood, mcmc.getOperatorSchedule(), logger);

            if (!xo.getAttribute(SPAWN,true))
                mle.setSpawnable(false);

            return mle;
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

        private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
                AttributeRule.newIntegerRule(CHAIN_LENGTH),
                AttributeRule.newIntegerRule(PATH_STEPS),
                AttributeRule.newBooleanRule(SPAWN,true),
                new ElementRule(MCMC.class),
                new ElementRule(PathLikelihood.class),
                new ElementRule(MCLogger.class)
        };

    };

    public String getId() {
        return id;
    }

    public void setId(String id) {
        this.id = id;
    }

    // PRIVATE TRANSIENTS

    /**
     * this markov chain does most of the work.
     */
    private MarkovChain mc;

    private String id = null;

    private int currentState;
    private int chainLength;

    private int burnin;
    private int pathSteps;
    private double pathDelta;
    private double pathParameter;

    private MCLogger logger;

    private PathLikelihood pathLikelihood;

    public static final String MARGINAL_LIKELIHOOD_ESTIMATOR = "marginalLikelihoodEstimator";
    public static final String CHAIN_LENGTH = "chainLength";
    public static final String PATH_STEPS = "pathSteps";
    public static final String SPAWN ="spawn";
}