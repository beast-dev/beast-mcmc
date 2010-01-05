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
import dr.inference.operators.CombinedOperatorSchedule;
import dr.inference.operators.OperatorSchedule;
import dr.inference.prior.Prior;
import dr.util.Identifiable;
import dr.xml.*;

/**
 * An MCMC analysis that estimates parameters of a probabilistic model.
 *
 * @author Andrew Rambaut
 * @author Alex Alekseyenko
 * @version $Id: MCMC.java,v 1.41 2005/07/11 14:06:25 rambaut Exp $
 */
public class MarginalLikelihoodEstimator implements Runnable, Identifiable {

    public MarginalLikelihoodEstimator(String id, int chainLength, int burninLength, int pathSteps,
                                       boolean linear, boolean lacing,
                                       PathLikelihood pathLikelihood,
                                       OperatorSchedule schedule,
                                       MCLogger logger) {

        this.id = id;
        this.chainLength = chainLength;
        this.pathSteps = pathSteps;
        this.linear = linear;
        this.lacing = lacing;
        this.burninLength = burninLength;

        MCMCCriterion criterion = new MCMCCriterion();

        pathDelta = 1.0 / pathSteps;
        pathParameter = 1.0;

        this.pathLikelihood = pathLikelihood;
        pathLikelihood.setPathParameter(pathParameter);

        mc = new MarkovChain(Prior.UNIFORM_PRIOR, pathLikelihood, schedule, criterion, 0, 0, false);

        this.logger = logger;
    }

    public void linearIntegration() {
        if (burninLength == -1)
            burnin = (int) (0.1 * chainLength);
        else
            burnin = burninLength;
        mc.setCurrentLength(0);
        for (int step = 0; step < pathSteps; step++) {
            pathLikelihood.setPathParameter(pathParameter);
            reportIteration(pathParameter, chainLength, burnin);
            mc.runChain(chainLength + burnin, false/*, 0*/);
            pathParameter -= pathDelta;

        }

        pathLikelihood.setPathParameter(0.0);
        reportIteration(pathParameter, chainLength, burnin);
        mc.runChain(chainLength + burnin, false/*, 0*/);
    }

    public void linearShoeLacing() {
        if (burninLength == -1) {
            burnin = (int) (0.1 * chainLength);
        } else {
            burnin = burninLength;
        }
        mc.setCurrentLength(0);
        if (pathSteps % 2 == 0) //Works only for odd number of steps
            pathSteps += 1;
        for (int step = 0; step < pathSteps - 1; step++) {
            pathLikelihood.setPathParameter(pathParameter);
            reportIteration(pathParameter, chainLength, burnin);
            mc.runChain(chainLength + burnin, false/*, 0*/);
            if (step == 0)
                pathParameter -= 2 * pathDelta;
            else if (step % 2 == 0) {
                pathParameter -= 3 * pathDelta;
            } else {
                pathParameter += pathDelta;
            }

        }

        pathLikelihood.setPathParameter(0.0);
        reportIteration(pathParameter, chainLength, burnin);
        mc.runChain(chainLength + burnin, false/*, 0*/);

    }

    public void geometricShoeLacing() {
        //int logCount = chainLength / logger.getLogEvery();
        mc.setCurrentLength(0);
        if (pathSteps % 2 == 0) //Works only for odd number of steps
            pathSteps += 1;

        for (int step = 0; step < pathSteps; step++) {
            double p = 1.0 - (((double) step) / pathSteps);
            int cl = (int) (20.0 * chainLength / ((19.0 * p) + 1));
//            logger.setLogEvery(cl / logCount);
            if (burninLength == -1)
                burnin = (int) (0.1 * cl);
            else
                burnin = burninLength;

            pathLikelihood.setPathParameter(pathParameter);
            reportIteration(pathParameter, cl, burnin);
            mc.runChain(cl + burnin, false/*, 0*/);
            if (step == 0) {
                pathParameter /= 4;
            } else if (step % 2 == 0) {
                pathParameter /= 8;
            } else {
                pathParameter *= 2;
            }
        }

        int cl = 20 * chainLength;
//        logger.setLogEvery(cl / logCount);
        if (burninLength == -1)
            burnin = (int) (0.1 * cl);
        else
            burnin = burninLength;
        pathLikelihood.setPathParameter(0.0);
        reportIteration(pathParameter, cl, burnin);
        mc.runChain(cl + burnin, false/*, 0*/);
    }

    private void reportIteration(double pathParameter, int cl, int burn) {
        System.out.println("Attempting theta = " + pathParameter + " for " + cl + " iterations + " + burn + " burnin.");
    }

    public void geometricIntegration() {
        //int logCount = chainLength / logger.getLogEvery();
        mc.setCurrentLength(0);
        for (int step = 0; step < pathSteps; step++) {
            double p = 1.0 - (((double) step) / pathSteps);
            int cl = (int) (20.0 * chainLength / ((19.0 * p) + 1));
//            logger.setLogEvery(cl / logCount);
            if (burninLength == -1)
                burnin = (int) (0.1 * cl);
            else burnin = burninLength;

            pathLikelihood.setPathParameter(pathParameter);
            reportIteration(pathParameter, cl, burnin);
            mc.runChain(cl + burnin, false/*, 0*/);
            pathParameter /= 2;

        }

        int cl = 20 * chainLength;
//        logger.setLogEvery(cl / logCount);
        if (burninLength == -1)
            burnin = (int) (0.1 * cl);
        else
            burnin = burninLength;
        pathLikelihood.setPathParameter(0.0);
        reportIteration(pathParameter, cl, burnin);
        mc.runChain(cl + burnin, false/*, 0*/);
    }

    public void run() {

        logger.startLogging();
        mc.addMarkovChainListener(chainListener);

        if (linear) {
            if (lacing)
                linearShoeLacing();
            else
                linearIntegration();
        } else {
            if (lacing)
                geometricShoeLacing();
            else
                geometricIntegration();
        }

        mc.removeMarkovChainListener(chainListener);
    }

    private final MarkovChainListener chainListener = new MarkovChainListener() {

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

            PathLikelihood pathLikelihood = (PathLikelihood) xo.getChild(PathLikelihood.class);
            MCLogger logger = (MCLogger) xo.getChild(MCLogger.class);

            int chainLength = xo.getIntegerAttribute(CHAIN_LENGTH);
            int pathSteps = xo.getIntegerAttribute(PATH_STEPS);
            int burninLength = -1;
            if (xo.hasAttribute(BURNIN)) {
                burninLength = xo.getIntegerAttribute(BURNIN);
            }
            boolean linear = true;
            boolean lacing = false;

            if (!xo.getAttribute(LINEAR, true))
                linear = false;

            if (xo.getAttribute(LACING, true))
                lacing = true;


            for (int i = 0; i < xo.getChildCount(); i++) {
                Object child = xo.getChild(i);
                if (child instanceof Logger) {
                }
            }

            CombinedOperatorSchedule os = new CombinedOperatorSchedule();

            XMLObject mcmcXML = xo.getChild(MCMC);
            for (int i = 0; i < mcmcXML.getChildCount(); ++i) {
                if (mcmcXML.getChild(i) instanceof MCMC) {
                    MCMC mcmc = (MCMC) mcmcXML.getChild(i);
                    os.addOperatorSchedule(mcmc.getOperatorSchedule());
                }
            }

            if (os.getScheduleCount() == 0) {
                System.err.println("Error: no mcmc objects provided in construction. Bayes Factor estimation will likely fail.");
            }

            java.util.logging.Logger.getLogger("dr.inference").info("Creating the Marginal Likelihood Estimator chain:" +
                    "\n  chainLength=" + chainLength +
                    "\n  pathSteps=" + pathSteps);

            MarginalLikelihoodEstimator mle = new MarginalLikelihoodEstimator(MARGINAL_LIKELIHOOD_ESTIMATOR, chainLength,
                    burninLength, pathSteps, linear, lacing, pathLikelihood, os, logger);

            if (!xo.getAttribute(SPAWN, true))
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

        private final XMLSyntaxRule[] rules = {
                AttributeRule.newIntegerRule(CHAIN_LENGTH),
                AttributeRule.newIntegerRule(PATH_STEPS),
                AttributeRule.newIntegerRule(BURNIN, true),
                AttributeRule.newBooleanRule(LINEAR, true),
                AttributeRule.newBooleanRule(LACING, true),
                AttributeRule.newBooleanRule(SPAWN, true),
                new ElementRule(MCMC,
                        new XMLSyntaxRule[]{new ElementRule(MCMC.class, 1, Integer.MAX_VALUE)}, false),
                //new ElementRule(MCMC.class),
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
    private final MarkovChain mc;

    private String id = null;

    private int currentState;
    private final int chainLength;

    private int burnin;
    private final int burninLength;
    private int pathSteps;
    private final boolean linear;
    private final boolean lacing;
    private final double pathDelta;
    private double pathParameter;

    private final MCLogger logger;

    private final PathLikelihood pathLikelihood;

    public static final String MARGINAL_LIKELIHOOD_ESTIMATOR = "marginalLikelihoodEstimator";
    public static final String CHAIN_LENGTH = "chainLength";
    public static final String PATH_STEPS = "pathSteps";
    public static final String LINEAR = "linear";
    public static final String LACING = "lacing";
    public static final String SPAWN = "spawn";
    public static final String BURNIN = "burnin";
    public static final String MCMC = "samplers";
}