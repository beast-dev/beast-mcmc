/*
 * SkylineLikelihood.java
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

package dr.evomodel.coalescent;

import dr.evolution.coalescent.ConstantPopulation;
import dr.evolution.tree.Tree;
import dr.evolution.util.Units;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.Likelihood;
import dr.inference.model.Parameter;
import dr.xml.*;
import no.uib.cipr.matrix.SymmTridiagMatrix;

/**
 * A likelihood function for the Gaussian Markov random field population trajectory.
 * *
 *
 * @author Erik Bloomquist
 * @author Vladimir Minin
 * @author Marc Suchard
 * @version $Id: GMRFSkylineLikelihood.java,v 1.3 2007/03/20 22:40:04 msuchard Exp $
 */
public class GMRFSkylineLikelihood extends CoalescentLikelihood {

    // PUBLIC STUFF

    public static final String SKYLINE_LIKELIHOOD = "gmrfSkyLineLikelihood";
    public static final String POPULATION_PARAMETER = "populationSizes";
    public static final String PRECISION_PARAMETER = "precisionParameter";
    public static final String POPULATION_TREE = "populationTree";
    public static final String LAMBDA_PARAMETER = "lambdaParameter";

    // PRIVATE STUFF

    protected Parameter popSizeParameter;
    protected Parameter precisionParameter;
    protected Parameter lambdaParameter;
    protected double[] gmrfWeights;
    protected int fieldLength;
    protected double[] coalescentIntervals;
    protected double[] storedCoalescentIntervals;
    protected double[] sufficientStatistics;
    protected double[] storedSufficientStatistics;

    protected SymmTridiagMatrix weightMatrix;
    protected SymmTridiagMatrix storedWeightMatrix;

    public GMRFSkylineLikelihood() {
        super(SKYLINE_LIKELIHOOD);
    }

    public GMRFSkylineLikelihood(Tree tree, Parameter popParameter, Parameter precParameter, Parameter lambda) {
        super(SKYLINE_LIKELIHOOD);

        this.popSizeParameter = popParameter;
        this.precisionParameter = precParameter;
        this.lambdaParameter = lambda;
        int tips = tree.getExternalNodeCount();
        fieldLength = popSizeParameter.getDimension();
        if (tips - fieldLength != 1) {
            throw new IllegalArgumentException("Number of tips (" + tips + ") must be one greater than number of pop sizes (" + fieldLength + ")");
        }

        this.tree = tree;
        if (tree instanceof TreeModel) {
            addModel((TreeModel) tree);
        }
        addParameter(popSizeParameter);
        addParameter(precisionParameter);
        addParameter(lambdaParameter);
        setupIntervals();
        coalescentIntervals = new double[fieldLength];
        storedCoalescentIntervals = new double[fieldLength];
        sufficientStatistics = new double[fieldLength];
        storedSufficientStatistics = new double[fieldLength];

        setupGMRFWeights();

        addStatistic(new DeltaStatistic());
    }

    // **************************************************************
    // Likelihood IMPLEMENTATION
    // **************************************************************

    public double getLogLikelihood() {
        if (!likelihoodKnown) {
            logLikelihood = calculateLogLikelihood();
            likelihoodKnown = true;
        }
        return logLikelihood;
    }

    public double[] getSufficientStatistics() {
        return sufficientStatistics;
    }

    protected void setupGMRFWeights() {

        int index = 0;

        double length = 0;
        double weight = 0;
        for (int i = 0; i < getIntervalCount(); i++) {
            length += getInterval(i);
            weight += getInterval(i) * getLineageCount(i) * (getLineageCount(i) - 1);
            if (getIntervalType(i) == CoalescentEventType.COALESCENT) {
                coalescentIntervals[index] = length;
                sufficientStatistics[index] = weight / 2.0;
                index++;
                length = 0;
                weight = 0;
            }
        }

        //Set up the weight Matrix
        double[] offdiag = new double[fieldLength - 1];
        double[] diag = new double[fieldLength];

//        double precision = precisionParameter.getParameterValue(0);

        //First set up the offdiagonal entries;
        for (int i = 0; i < fieldLength - 1; i++) {
            offdiag[i] = -2.0 / (coalescentIntervals[i] + coalescentIntervals[i + 1]);
        }

        //Then set up the diagonal entries;
        for (int i = 1; i < fieldLength - 1; i++)
            diag[i] = -(offdiag[i] + offdiag[i - 1]);

        //Take care of the endpoints
        diag[0] = -offdiag[0];
        diag[fieldLength - 1] = -offdiag[fieldLength - 2];


        weightMatrix = new SymmTridiagMatrix(diag, offdiag);

    }


    public SymmTridiagMatrix getScaledWeightMatrix(double precision) {
        SymmTridiagMatrix a = weightMatrix.copy();
        for (int i = 0; i < a.numRows() - 1; i++) {
            a.set(i, i, a.get(i, i) * precision);
            a.set(i + 1, i, a.get(i + 1, i) * precision);
        }
        a.set(fieldLength - 1, fieldLength - 1, a.get(fieldLength - 1, fieldLength - 1) * precision);
        return a;
    }

    public SymmTridiagMatrix getCopyWeightMatrix() {
        return weightMatrix.copy();
    }


    protected void storeState() {
        super.storeState();
        System.arraycopy(coalescentIntervals, 0, storedCoalescentIntervals, 0, coalescentIntervals.length);
        System.arraycopy(sufficientStatistics, 0, storedSufficientStatistics, 0, sufficientStatistics.length);
        storedWeightMatrix = weightMatrix.copy();
    }


    protected void restoreState() {
        super.restoreState();
        System.arraycopy(storedCoalescentIntervals, 0, coalescentIntervals, 0, storedCoalescentIntervals.length);
        System.arraycopy(storedSufficientStatistics, 0, sufficientStatistics, 0, storedSufficientStatistics.length);
        weightMatrix = storedWeightMatrix;

    }


    protected void handleParameterChangedEvent(Parameter parameter, int index) {
        likelihoodKnown = false;
        // Parameters (precision and popsizes do not change intervals or GMRF sufficient statistics

        // todo precision and lambda both change Q

    }

    /**
     * Calculates the log likelihood of this set of coalescent intervals,
     * given a demographic model.
     */
    public double calculateLogLikelihood() {

        if (!intervalsKnown) {
            // intervalsKnown -> false when handleModelChanged event occurs in super.
            setupIntervals();
            setupGMRFWeights();
        }

        double logL = 0.0;

        double currentTime = 0.0;

        int popIndex = 0;

        ConstantPopulation cp = new ConstantPopulation(Units.Type.YEARS);

        for (int j = 0; j < intervalCount; j++) {

            cp.setN0(Math.exp(popSizeParameter.getParameterValue(popIndex)));
            if (getIntervalType(j) == CoalescentEventType.COALESCENT) {
                popIndex += 1;
            }

            logL += calculateIntervalLikelihood(cp, intervals[j], currentTime, lineageCounts[j], getIntervalType(j));

            // insert zero-length coalescent intervals
//			int diff = getCoalescentEvents(j) - 1;
//			System.err.println("Diff = "+diff);
//			for (int k = 0; k < diff; k++) {
//				cp.setN0(popSizeParameter.getParameterValue(popIndex));
//				logL += calculateIntervalLikelihood(cp, 0.0, currentTime, lineageCounts[j] - k - 1, CoalescentEventType.COALESCENT);
//				popIndex += 1;
//			}

            currentTime += intervals[j];

        }

        // Calculate GMRF density; here GMRF = RW(1)
        logL += -0.5 * calculateWeightedSSE() * precisionParameter.getParameterValue(0);

/*

		WinBUGS code to fixed tree:  (A:4.0,(B:2.0,(C:0.5,D:1.0):1.0):2.0)

		model {

			stat1 ~ dexp(rate[1])
			stat2 ~ dexp(rate[2])
			stat3 ~ dexp(rate[3])

			rate[1] <- 1 / exp(theta[1])
			rate[2] <- 1 / exp(theta[2])
			rate[3] <- 1 / exp(theta[3])

			theta[1] ~ dnorm(0, 0.001)
			theta[2] ~ dnorm(theta[1], weight[1])
			theta[3] ~ dnorm(theta[2], weight[2])

			weight[1] <- tau / 1.0
			weight[2] <- tau / 1.5

			tau ~ dgamma(1,0.3333)

			stat1 <- 9 / 2
			stat2 <- 6 / 2
			stat3 <- 4 / 2

		} 

*/

//		return logL;
        return 0;
    }


    public Parameter getPrecisionParameter() {
        return precisionParameter;
    }

    public Parameter getPopSizeParameter() {
        return popSizeParameter;
    }

    public Parameter getLambdaParameter() {
        return lambdaParameter;
    }

    public SymmTridiagMatrix getWeightMatrix() {
        return weightMatrix.copy();
    }


    public double calculateWeightedSSE() {
        double weightedSSE = 0;
        double currentPopSize = popSizeParameter.getParameterValue(0);
        double currentInterval = coalescentIntervals[0];
        for (int j = 1; j < fieldLength; j++) {
            double nextPopSize = popSizeParameter.getParameterValue(j);
            double nextInterval = coalescentIntervals[j];
            double delta = nextPopSize - currentPopSize;
            double weight = (currentInterval + nextInterval) / 2.0;
            weightedSSE += delta * delta / weight;
            currentPopSize = nextPopSize;
            currentInterval = nextInterval;
        }
        return weightedSSE;

    }

    // ****************************************************************
    // Private and protected stuff
    // ****************************************************************

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return SKYLINE_LIKELIHOOD;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            XMLObject cxo = (XMLObject) xo.getChild(POPULATION_PARAMETER);
            Parameter popParameter = (Parameter) cxo.getChild(Parameter.class);

            cxo = (XMLObject) xo.getChild(PRECISION_PARAMETER);
            Parameter precParameter = (Parameter) cxo.getChild(Parameter.class);

            cxo = (XMLObject) xo.getChild(POPULATION_TREE);
            TreeModel treeModel = (TreeModel) cxo.getChild(TreeModel.class);

            cxo = (XMLObject) xo.getChild(LAMBDA_PARAMETER);
            Parameter lambda = (Parameter) cxo.getChild(Parameter.class);

            return new GMRFSkylineLikelihood(treeModel, popParameter, precParameter, lambda);
        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public String getParserDescription() {
            return "This element represents the likelihood of the tree given the population size vector.";
        }

        public Class getReturnType() {
            return Likelihood.class;
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
                new ElementRule(POPULATION_PARAMETER, new XMLSyntaxRule[]{
                        new ElementRule(Parameter.class)
                }),
                new ElementRule(PRECISION_PARAMETER, new XMLSyntaxRule[]{
                        new ElementRule(Parameter.class)
                }),
                new ElementRule(POPULATION_TREE, new XMLSyntaxRule[]{
                        new ElementRule(TreeModel.class)
                }),
                new ElementRule(LAMBDA_PARAMETER, new XMLSyntaxRule[]{
                        new ElementRule(Parameter.class)
                }),
        };
    };


}