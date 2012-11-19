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
import dr.inference.operators.OperatorAnalysisPrinter;
import dr.inference.operators.OperatorSchedule;
import dr.inference.prior.Prior;
import dr.util.Identifiable;
import dr.xml.*;
import org.apache.commons.math.MathException;
import org.apache.commons.math.distribution.BetaDistributionImpl;

/**
 * An MCMC analysis that estimates parameters of a probabilistic model.
 *
 * @author Andrew Rambaut
 * @author Alex Alekseyenko
 * @version $Id: MCMC.java,v 1.41 2005/07/11 14:06:25 rambaut Exp $
 */
public class MarginalLikelihoodEstimator implements Runnable, Identifiable {

    public MarginalLikelihoodEstimator(String id, int chainLength, int burninLength, int pathSteps, double fixedRunValue,
//                                       boolean linear, boolean lacing,
PathScheme scheme,
PathLikelihood pathLikelihood,
OperatorSchedule schedule,
MCLogger logger) {

        this.id = id;
        this.chainLength = chainLength;
        this.pathSteps = pathSteps;
        this.scheme = scheme;
        this.schedule = schedule;
        this.fixedRunValue = fixedRunValue;
        // deprecated
        // this.linear = (scheme == PathScheme.LINEAR);
        // this.lacing = false; // Was not such a good idea

        this.burninLength = burninLength;

        MCMCCriterion criterion = new MCMCCriterion();

        pathDelta = 1.0 / pathSteps;
        pathParameter = 1.0;

        this.pathLikelihood = pathLikelihood;
        pathLikelihood.setPathParameter(pathParameter);

        mc = new MarkovChain(Prior.UNIFORM_PRIOR, pathLikelihood, schedule, criterion, 0, 0, true);

        this.logger = logger;
    }

    private void setDefaultBurnin() {
        if (burninLength == -1) {
            burnin = (int) (0.1 * chainLength);
        } else {
            burnin = burninLength;
        }
    }

    public void integrate(Integrator scheme) {
        setDefaultBurnin();
        mc.setCurrentLength(burnin);
        scheme.init();
        for (pathParameter = scheme.nextPathParameter(); pathParameter >= 0; pathParameter = scheme.nextPathParameter()) {
            pathLikelihood.setPathParameter(pathParameter);
            reportIteration(pathParameter, chainLength, burnin, scheme.pathSteps, scheme.step);
            long cl = mc.getCurrentLength();
            mc.setCurrentLength(0);
            mc.runChain(burnin, false/*, 0*/);
            mc.setCurrentLength(cl);
            mc.runChain(chainLength, false);
            (new OperatorAnalysisPrinter(schedule)).showOperatorAnalysis(System.out);
            ((CombinedOperatorSchedule) schedule).reset();
        }
    }

    public abstract class Integrator {
        protected int step;
        protected int pathSteps;

        protected Integrator(int pathSteps) {
            this.pathSteps = pathSteps;
        }

        public void init() {
            step = 0;
        }

        abstract double nextPathParameter();
    }
    
    public class FixedThetaRun extends Integrator {
    	private double value;
    	
    	public FixedThetaRun(double value) {
    		super(1);
    		this.value = value;
    	}
    	
    	double nextPathParameter() {
    		if (step == 0) {
    			step++;
    			return value;
    		} else {
    			return -1.0;
    		}
    	}
    	
    }

    public class LinearIntegrator extends Integrator {
        public LinearIntegrator(int pathSteps) {
            super(pathSteps);
        }

        double nextPathParameter() {
            if (step > pathSteps) {
                return -1;
            }
            double pathParameter = 1.0 - (double)step / (double)(pathSteps - 1);
            step = step + 1;
            return pathParameter;
        }
    }

    public class SigmoidIntegrator extends Integrator {
    	private double alpha;

    	public SigmoidIntegrator(double alpha, int pathSteps) {
    		super(pathSteps);
    		this.alpha = alpha;
    	}

    	double nextPathParameter() {
    		if (step == 0) {
    			step++;
    			return 1.0;
    		} else if (step == pathSteps) {
    			step++;
    			return 0.0;
    		} else if (step > pathSteps) {
    			return -1.0;
    		} else {
    			double xvalue = ((pathSteps - step)/((double)pathSteps)) - 0.5;
    			step++;
    			return Math.exp(alpha*xvalue)/(Math.exp(alpha*xvalue) + Math.exp(-alpha*xvalue));
    		}
    	}
    }

    public class BetaQuantileIntegrator extends Integrator {
    	private double alpha;

    	public BetaQuantileIntegrator(double alpha, int pathSteps) {
    		super(pathSteps);
    		this.alpha = alpha;
    	}

    	double nextPathParameter() {
    		if (step > pathSteps)
                return -1;
    		double result = Math.pow((pathSteps - step)/((double)pathSteps), 1.0/alpha);
    		step++;
    		return result;
    	}
    }

    public class BetaIntegrator extends Integrator {
        private BetaDistributionImpl betaDistribution;

        public BetaIntegrator(double alpha, double beta, int pathSteps) {
            super(pathSteps);
            this.betaDistribution = new BetaDistributionImpl(alpha, beta);
        }

        double nextPathParameter() {
            if (step > pathSteps)
                return -1;
            if (step == 0) {
                step += 1;
                return 1.0;
            } else if (step + 1 < pathSteps) {
                double ratio = (double) step / (double) (pathSteps - 1);
                try {
                    step += 1;
                    return 1.0 - betaDistribution.inverseCumulativeProbability(ratio);
                } catch (MathException e) {
                    e.printStackTrace();
                }
            }
            step += 1;
            return 0.0;
        }
    }

    public class GeometricIntegrator extends Integrator {

        public GeometricIntegrator(int pathSteps) {
            super(pathSteps);
        }

        double nextPathParameter() {
            if (step > pathSteps) {
                return -1;
            }
            if (step == pathSteps) { //pathSteps instead of pathSteps - 1
                step += 1;
                return 0;
            }

            step += 1;
            return Math.pow(2, -(step - 1));
        }
    }

    /*public void linearIntegration() {
        setDefaultBurnin();
        mc.setCurrentLength(0);
        for (int step = 0; step < pathSteps; step++) {
            pathLikelihood.setPathParameter(pathParameter);
            reportIteration(pathParameter, chainLength, burnin);
            //mc.runChain(chainLength + burnin, false, 0);
            mc.runChain(chainLength + burnin, false);
            pathParameter -= pathDelta;
        }
        pathLikelihood.setPathParameter(0.0);
        reportIteration(pathParameter, chainLength, burnin);
        //mc.runChain(chainLength + burnin, false, 0);
        mc.runChain(chainLength + burnin, false);
    }*/

    /*public void betaIntegration(double alpha, double beta) {
        setDefaultBurnin();
        mc.setCurrentLength(0);

        BetaDistributionImpl betaDistribution = new BetaDistributionImpl(alpha, beta);

        for (int step = 0; step < pathSteps; step++) {
            if (step == 0) {
                pathParameter = 1.0;
            } else if (step + 1 < pathSteps) {
                double ratio = (double) step / (double) (pathSteps - 1);
                try {
                    pathParameter = 1.0 - betaDistribution.inverseCumulativeProbability(ratio);
                } catch (MathException e) {
                    e.printStackTrace();
                }
            } else {
                pathParameter = 0.0;
            }
            pathLikelihood.setPathParameter(pathParameter);
            reportIteration(pathParameter, chainLength, burnin);
            //mc.runChain(chainLength + burnin, false, 0);
            mc.runChain(chainLength + burnin, false);
            (new OperatorAnalysisPrinter(schedule)).showOperatorAnalysis(System.out);
            ((CombinedOperatorSchedule) schedule).reset();
        }
    }*/

    private void reportIteration(double pathParameter, long chainLength, long burnin, long totalSteps, long steps) {
        System.out.println("Attempting theta ("+steps+"/" + totalSteps +") = " + pathParameter + " for " + chainLength + " iterations + " + burnin + " burnin.");
    }

    public void run() {

        logger.startLogging();
        mc.addMarkovChainListener(chainListener);

        /*switch (scheme) {
            case LINEAR:
                linearIntegration();
                break;
            case OLD_GEOMETRIC:
                geometricIntegration();
                break;
            case ONE_SIDED_BETA:
                betaIntegration(1.0, betaFactor);
                break;
            case BETA:
                betaIntegration(alphaFactor, betaFactor);
                break;
            default:
                throw new RuntimeException("Illegal path scheme");
        }*/

        switch (scheme) {
        	case FIXED:
        		integrate(new FixedThetaRun(fixedRunValue));
        		break;
            case LINEAR:
                integrate(new LinearIntegrator(pathSteps));
                break;
            case GEOMETRIC:
                integrate(new GeometricIntegrator(pathSteps));
                break;
            case ONE_SIDED_BETA:
                integrate(new BetaIntegrator(1.0, betaFactor, pathSteps));
                break;
            case BETA:
                integrate(new BetaIntegrator(alphaFactor, betaFactor, pathSteps));
                break;
            case BETA_QUANTILE:
            	integrate(new BetaQuantileIntegrator(alphaFactor, pathSteps));
            	break;
            case SIGMOID:
            	integrate(new SigmoidIntegrator(alphaFactor, pathSteps));
            	break;
            default:
                throw new RuntimeException("Illegal path scheme");
        }

        mc.removeMarkovChainListener(chainListener);
    }

    private final MarkovChainListener chainListener = new MarkovChainListener() {

        // MarkovChainListener interface *******************************************
        // for receiving messages from subordinate MarkovChain

        /**
         * Called to update the current model keepEvery states.
         */
        public void currentState(long state, Model currentModel) {

            currentState = state;

            if (currentState >= burnin) {
                logger.log(state);
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
            (new OperatorAnalysisPrinter(schedule)).showOperatorAnalysis(System.out);
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

    public void setAlphaFactor(double alpha) {
        alphaFactor = alpha;
    }

    public void setBetaFactor(double beta) {
        betaFactor = beta;
    }

    public double getAlphaFactor() {
        return alphaFactor;
    }

    public double getBetaFactor() {
        return betaFactor;
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

            int prerunLength = -1;
            if (xo.hasAttribute(PRERUN)) {
                prerunLength = xo.getIntegerAttribute(PRERUN);
            }
            double fixedRunValue = -1.0;
            if (xo.hasAttribute(FIXED_VALUE)) {
            	fixedRunValue = xo.getDoubleAttribute(FIXED_VALUE);
            }

            // deprecated
            boolean linear = xo.getAttribute(LINEAR, true);
            // boolean lacing = xo.getAttribute(LACING,false);
            PathScheme scheme;
            if (linear) {
                scheme = PathScheme.LINEAR;
            } else {
                scheme = PathScheme.GEOMETRIC;
            }

            // new approach
            if (xo.hasAttribute(PATH_SCHEME)) { // change to: getAttribute once deprecated approach removed
                scheme = PathScheme.parseFromString(xo.getAttribute(PATH_SCHEME, PathScheme.LINEAR.getText()));
            }

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
                    if (prerunLength > 0) {
                        java.util.logging.Logger.getLogger("dr.inference").info("Path Sampling Marginal Likelihood Estimator:\n\tEquilibrating chain " + mcmc.getId() + " for " + prerunLength + " iterations.");
                        for (Logger log : mcmc.getLoggers()) { // Stop the loggers, so nothing gets written to normal output
                            log.stopLogging();
                        }
                        mcmc.getMarkovChain().runChain(prerunLength, false);
                    }
                    os.addOperatorSchedule(mcmc.getOperatorSchedule());
                }
            }

            if (os.getScheduleCount() == 0) {
                System.err.println("Error: no mcmc objects provided in construction. Bayes Factor estimation will likely fail.");
            }

            MarginalLikelihoodEstimator mle = new MarginalLikelihoodEstimator(MARGINAL_LIKELIHOOD_ESTIMATOR, chainLength,
                    burninLength, pathSteps, fixedRunValue, scheme, pathLikelihood, os, logger);

            if (!xo.getAttribute(SPAWN, true))
                mle.setSpawnable(false);

            if (xo.hasAttribute(ALPHA)) {
                mle.setAlphaFactor(xo.getAttribute(ALPHA, 0.5));
            }

            if (xo.hasAttribute(BETA)) {
                mle.setBetaFactor(xo.getAttribute(BETA, 0.5));
            }

            String alphaBetaText = "";
            if (scheme == PathScheme.ONE_SIDED_BETA) {
                alphaBetaText += "(1," + mle.getBetaFactor() + ")";
            } else if (scheme == PathScheme.BETA) {
                alphaBetaText += "(" + mle.getAlphaFactor() + "," + mle.getBetaFactor() + ")";
            } else if (scheme == PathScheme.BETA_QUANTILE) {
            	alphaBetaText += "(" + mle.getAlphaFactor() + ")";
            } else if (scheme == PathScheme.SIGMOID) {
            	alphaBetaText += "(" + mle.getAlphaFactor() + ")";
            }
            java.util.logging.Logger.getLogger("dr.inference").info("\nCreating the Marginal Likelihood Estimator chain:" +
                    "\n  chainLength=" + chainLength +
                    "\n  pathSteps=" + pathSteps +
                    "\n  pathScheme=" + scheme.getText() + alphaBetaText +
                    "\n  If you use these results, please cite:" +
                    "\n    Guy Baele, Philippe Lemey, Trevor Bedford, Andrew Rambaut, Marc A. Suchard, and Alexander V. Alekseyenko." +
                    "\n    2012. Improving the accuracy of demographic and molecular clock model comparison while accommodating " + 
                    "\n          phylogenetic uncertainty. Mol. Biol. Evol. (in press).");
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
                AttributeRule.newIntegerRule(PRERUN, true),
                AttributeRule.newBooleanRule(LINEAR, true),
                AttributeRule.newBooleanRule(LACING, true),
                AttributeRule.newBooleanRule(SPAWN, true),
                AttributeRule.newStringRule(PATH_SCHEME, true),
                AttributeRule.newDoubleRule(FIXED_VALUE, true),
                AttributeRule.newDoubleRule(ALPHA, true),
                AttributeRule.newDoubleRule(BETA, true),
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

    enum PathScheme {
    	FIXED("fixed"),
        LINEAR("linear"),
        GEOMETRIC("geometric"),
        BETA("beta"),
        ONE_SIDED_BETA("oneSidedBeta"),
        BETA_QUANTILE("betaQuantile"),
        SIGMOID("sigmoid");

        PathScheme(String text) {
            this.text = text;
        }

        public String getText() {
            return text;
        }

        private final String text;

        public static PathScheme parseFromString(String text) {
            for (PathScheme scheme : PathScheme.values()) {
                if (scheme.getText().compareToIgnoreCase(text) == 0)
                    return scheme;
            }
            return null;
        }
    }


    // PRIVATE TRANSIENTS

    /**
     * this markov chain does most of the work.
     */
    private final MarkovChain mc;
    private OperatorSchedule schedule;

    private String id = null;

    private long currentState;
    private final long chainLength;

    private long burnin;
    private final long burninLength;
    private int pathSteps;
    //    private final boolean linear;
    //    private final boolean lacing;
    private final PathScheme scheme;
    private double alphaFactor = 0.5;
    private double betaFactor = 0.5;
    private double fixedRunValue = -1.0;
    private final double pathDelta;
    private double pathParameter;

    private final MCLogger logger;

    private final PathLikelihood pathLikelihood;

    public static final String MARGINAL_LIKELIHOOD_ESTIMATOR = "marginalLikelihoodEstimator";
    public static final String CHAIN_LENGTH = "chainLength";
    public static final String PATH_STEPS = "pathSteps";
    public static final String FIXED = "fixed";
    public static final String LINEAR = "linear";
    public static final String LACING = "lacing";
    public static final String SPAWN = "spawn";
    public static final String BURNIN = "burnin";
    public static final String MCMC = "samplers";
    public static final String PATH_SCHEME = "pathScheme";
    public static final String FIXED_VALUE = "fixedValue";
    public static final String ALPHA = "alpha";
    public static final String BETA = "beta";
    public static final String PRERUN = "prerun";
}