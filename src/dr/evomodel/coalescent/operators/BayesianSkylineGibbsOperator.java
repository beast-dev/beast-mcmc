/*
 * BayesianSkylineGibbsOperator.java
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

package dr.evomodel.coalescent.operators;

import dr.evolution.coalescent.ConstantPopulation;
import dr.evolution.util.Units;
import dr.evomodel.coalescent.BayesianSkylineLikelihood;
import dr.evomodel.coalescent.OldAbstractCoalescentLikelihood;
import dr.evomodelxml.coalescent.operators.BayesianSkylineGibbsOperatorParser;
import dr.inference.model.Parameter;
import dr.inference.operators.SimpleMCMCOperator;
import dr.math.distributions.ExponentialDistribution;
import dr.math.distributions.GammaDistribution;

/**
 * A Gibbs operator for the Bayesian Skyline model, inspired by the Hey &
 * Nielsen PNAS paper.
 * <p/>
 * The operator is implemented as a standard proposal generator, returning a
 * Hastings ratio. The advantage of this scheme is that additional priors can be
 * included without modifying the operator. This is only be efficient for 'mild'
 * priors, since acceptance rates easily get very low because this operator
 * changes many parameters simultaneously.
 * <p/>
 * The proposal generator can be modified by providing upper and lower bounds.
 * Four types of prior are implemented, so that the proposal behaves as a true
 * Gibbs operator if the chosen prior matches the actual prior.
 * <p/>
 * jeffreys= exponentialMarkov=
 * false     false 				Uniform prior
 * true 	 false 				Jeffrey's prior on the geometric average population size
 * false 	 true 				An exponential Markov prior (as ExponentialMarkov)
 * true      true 				Same, with a Jeffrey's prior on the most recent pop size
 * <p/>
 * In addition, if reverse="true", the exponential Markov prior works in reverse
 * (any Jeffrey's prior applies to the last pop size parameter, and preceding
 * pop sizes are exponentially distributed conditional on the successor)
 * <p/>
 * The 'shape' parameter determines the shape parameter of the Gamma distribution used in the "exponential"
 * Markov prior.  The default value is 1.0, corresponding to an exponential distribution; otherwise a Gamma
 * distribution with that shape parameter, and a mean equal to the previous population size, is used.
 * <p/>
 * The 'iterations' parameter determines how many single-population-size Gibbs moves
 * will be executed per operation.  This should be set to once or twice the number
 * of dimensions.
 * <p/>
 * TODO: It's easy to extend this to use Gamma Markov priors (either on Ne or on 1/Ne)
 * instead of an exponential prior
 *
 * @author Gerton Lunter
 */

public class BayesianSkylineGibbsOperator extends SimpleMCMCOperator {

    public static final int STEPWISE_TYPE = 0;
    public static final int LINEAR_TYPE = 1;
    public static final int EXPONENTIAL_TYPE = 2;

    private BayesianSkylineLikelihood bayesianSkylineLikelihood;
    private Parameter populationSizeParameter;
    private Parameter groupSizeParameter;

    private double upperBound;
    private double lowerBound;

    private boolean jeffreysGeometricPrior; /* Jeffrey's prior on the geometric average */
    private boolean jeffreysPrior; /* Jeffrey's prior on the first Ne */
    private boolean exponentialMarkovPrior; /* Exponential or Gamma Markov prior */
    private double shape; /* shape parameter of the Gamma prior (1.0 for exponential prior) */
    private boolean reverse; /*
								 * If true, priors are parameterized by
								 * ancestral rather than descendant Ne
								 */

    private int iterations; /* Number of Gibbs moves to make per operation */

    public double getTargetAcceptanceProbability() {
        return 1.00;
    }

    public double getMinimumAcceptanceLevel() {
        return 0.99;
    }

    public double getMaximumAcceptanceLevel() {
        return 1.01;
    }

    public double getMinimumGoodAcceptanceLevel() {
        return 0.99;
    }

    public double getMaximumGoodAcceptanceLevel() {
        return 1.01;
    }

    public BayesianSkylineGibbsOperator(
            BayesianSkylineLikelihood bayesianSkylineLikelihood,
            Parameter populationSizeParameter, Parameter groupSizeParameter,
            int type, double weight, double lowerBound, double upperBound,
            boolean jeffreysPrior, boolean exponentialMarkov, double shape, boolean reverse,
            int iterations) {

        this.bayesianSkylineLikelihood = bayesianSkylineLikelihood;
        this.populationSizeParameter = populationSizeParameter;
        this.groupSizeParameter = groupSizeParameter;
        setWeight(weight);
        this.lowerBound = lowerBound;
        this.upperBound = upperBound;
        this.exponentialMarkovPrior = exponentialMarkov;
        this.shape = shape;
        this.iterations = iterations;
        if (exponentialMarkov) {
            this.jeffreysPrior = jeffreysPrior;
            this.jeffreysGeometricPrior = false;
        } else {
            this.jeffreysGeometricPrior = jeffreysPrior;
            this.jeffreysPrior = false;
        }
        this.reverse = reverse;

        assert (populationSizeParameter != null);
        assert (groupSizeParameter != null);
        assert (populationSizeParameter.getDimension() == groupSizeParameter
                .getDimension());
        assert (iterations >= 1);
        assert (lowerBound < upperBound);

        if (type != STEPWISE_TYPE) {
            throw new IllegalArgumentException(
                    "Should only get stepwise type here - sorry.");
        }

        System.out.println("Using a Bayesian skyline Gibbs operator (lo="
                + lowerBound + ", hi=" + upperBound + ", Jeffreys="
                + jeffreysPrior + ", exponentialMarkov=" + exponentialMarkov
                + " (shape=" + shape + "), reverse=" + reverse + ", iterations=" + iterations + ")");
    }

    public final String getPerformanceSuggestion() {
        return "";
    }

    // Samples one new population size, from a Gamma distribution with the given
    // shape and rate parameters,
    // and, when required, taking account of the exponential Markov prior on
    // either side
    double getSample(double exponents[], double gammaRates[], double popSizes[], int index) {

        int numGroups = popSizes.length;
        int direction = reverse ? -1 : 1;
        int firstIdx = reverse ? numGroups - 1 : 0;
        int lastIdx = reverse ? 0 : numGroups - 1;

        double exponent = exponents[index];
        double rate = gammaRates[index];
        double bias = 0.0;

        // Include a Jeffrey's prior on the geometric average of population
        // sizes
        if (jeffreysGeometricPrior) {
            exponent += 1.0 / numGroups;
        }

        // Include a Jeffrey's prior on the first population size
        if (jeffreysPrior && index == firstIdx) {
            exponent += 1.0;
        }

        // For a parameter other than the first, include an exponential prior
        // (1/Nprev) exp(-N/Nprev)
        // (The normalization factor is ignored, since it does not depend on N)
        //
        // Modified: include a Gamma prior
        //   N^(shape-1) exp(-N/scale) / scale^shape Gamma(shape)  with  scale = Nprev/shape
        // All factors not depending on N are ignored
        if (exponentialMarkovPrior && index != firstIdx) {
            //bias = popSizes[index - direction];
            bias = popSizes[index - direction] / shape;
            exponent += 1.0 - shape;
        }

        // Take into account the dependent prior (1/N) exp(-Nnext/N)
        // This time, the normalization factor, as well as the exponential
        // factor, must be included
        //
        // Modified: include a Gamma prior
        //   Nnext^(shape-1) exp(-Nnext*shape/N) shape^shape / N^shape Gamma(shape)
        if (exponentialMarkovPrior && index != lastIdx) {
            //exponent += 1.0;
            //rate = rate + popSizes[index + direction];
            exponent += shape;
            rate = rate + popSizes[index + direction] * shape;
        }

        // Check arguments
        // Note: the requirement that the exponent be > 1.0 (i.e. shape parameter > 0.0)
        //       could be relaxed when there is nonzero bias, but this is not implemented
        //       in GammaDistribution.nextExpGamma

        //if (exponent <= 1.0) {
        //    throw new IllegalArgumentException("Group size must be >= 2");
        //}

        // now sample
        double sample;
        int iters = 0;
        do {

            // now sample from (1/N)^exponent * exp(-rate/N) dN; this is
            // equivalent to sampling
            // from x^(exponent-2) * exp(-rate x) dx with x=1/N, which is the
            // Gamma distribution
            // with shape parameter exponent-1.
            if (bias > 0) {
                sample = 1.0 / GammaDistribution.nextExpGamma(exponent - 1.0, 1.0 / rate, bias);
            } else {
                sample = 1.0 / GammaDistribution.nextGamma(exponent - 1.0, 1.0 / rate);
            }
            iters += 1;

        } while ((sample < lowerBound || sample > upperBound) && iters < 100);

        if (iters == 100) {
            // fail
            return Double.NEGATIVE_INFINITY;
        }

        // calculate partial log likelihoods of old and new sample
        double oldLogL = (-exponent * Math.log(popSizes[index])) - rate
                / popSizes[index];
        double newLogL = (-exponent * Math.log(sample)) - rate / sample;

        if (bias > 0.0) {
            oldLogL += ExponentialDistribution.logPdf(popSizes[index], 1.0 / bias);
            newLogL += ExponentialDistribution.logPdf(sample, 1.0 / bias);
        }

        // store new sample
        popSizes[index] = sample;

        // return hastings ratio
        return oldLogL - newLogL;

    }

    public double doOperation() {

        if (!bayesianSkylineLikelihood.getIntervalsKnown())
            bayesianSkylineLikelihood.setupIntervals();

        assert (populationSizeParameter != null);
        assert (groupSizeParameter != null);

        // Now enter a loop similar to the likelihood calculation, except that
        // we now just collect
        // waiting times, and exponents of the population size parameters, for
        // each group.

        int groupIndex = 0;
        int subIndex = 0;
        int[] groupSizes = bayesianSkylineLikelihood.getGroupSizes();
        double[] groupEnds = bayesianSkylineLikelihood.getGroupHeights();

        ConstantPopulation cp = new ConstantPopulation(Units.Type.YEARS);
        double exponent = 0.0;
        double gammaRate = 0.0;
        double currentTime = 0.0;

        double[] populationSizes = new double[groupSizes.length];
        double[] exponents = new double[groupSizes.length];
        double[] gammaRates = new double[groupSizes.length];

        // first collect the shape and rate parameters of the Gamma
        // distributions describing the posteriors on N
        for (int j = 0; j < bayesianSkylineLikelihood.getIntervalCount(); j++) {

            // set the population size to the size of the middle of the current
            // interval
            double curpopsize = bayesianSkylineLikelihood.getPopSize(
                    groupIndex, currentTime + (bayesianSkylineLikelihood.getInterval(j) / 2.0),
                    groupEnds);
            populationSizes[groupIndex] = curpopsize;

            cp.setN0(curpopsize);

            exponent += bayesianSkylineLikelihood.calculateIntervalShapeParameter(cp,
                    bayesianSkylineLikelihood.getInterval(j),
                    currentTime,
                    bayesianSkylineLikelihood.getLineageCount(j),
                    bayesianSkylineLikelihood.getIntervalType(j));
            gammaRate += bayesianSkylineLikelihood.calculateIntervalRateParameter(cp,
                    bayesianSkylineLikelihood.getInterval(j),
                    currentTime,
                    bayesianSkylineLikelihood.getLineageCount(j),
                    bayesianSkylineLikelihood.getIntervalType(j))
                    * curpopsize;

            if (bayesianSkylineLikelihood.getIntervalType(j) == OldAbstractCoalescentLikelihood.CoalescentEventType.COALESCENT) {
                subIndex += 1;
                if (subIndex >= groupSizes[groupIndex]) {

                    exponents[groupIndex] = exponent;
                    gammaRates[groupIndex] = gammaRate;

                    // Reset accumulators
                    exponent = 0.0;
                    gammaRate = 0.0;

                    // Next group
                    groupIndex += 1;
                    subIndex = 0;
                }
            }

            // insert zero-length coalescent intervals (for multifurcating trees
            int diff = bayesianSkylineLikelihood.getCoalescentEvents(j) - 1;
            for (int k = 0; k < diff; k++) {

                curpopsize = bayesianSkylineLikelihood.getPopSize(groupIndex,
                        currentTime + (bayesianSkylineLikelihood.getInterval(j) / 2.0),
                        groupEnds);
                populationSizes[groupIndex] = curpopsize;

                cp.setN0(curpopsize);

                exponent += bayesianSkylineLikelihood.calculateIntervalShapeParameter(cp,
                        bayesianSkylineLikelihood.getInterval(j),
                        currentTime,
                        bayesianSkylineLikelihood.getLineageCount(j),
                        bayesianSkylineLikelihood.getIntervalType(j));
                gammaRate += bayesianSkylineLikelihood.calculateIntervalRateParameter(cp,
                        bayesianSkylineLikelihood.getInterval(j),
                        currentTime,
                        bayesianSkylineLikelihood.getLineageCount(j),
                        bayesianSkylineLikelihood.getIntervalType(j))
                        * curpopsize;

                subIndex += 1;
                if (subIndex >= groupSizes[groupIndex]) {

                    exponents[groupIndex] = exponent;
                    gammaRates[groupIndex] = gammaRate;

                    // Reset accumulators
                    exponent = 0.0;
                    gammaRate = 0.0;

                    // Next group
                    groupIndex += 1;
                    subIndex = 0;

                }
            }

            currentTime += bayesianSkylineLikelihood.getInterval(j);
        }

        // Next, sample new population sizes
        double hastingsRatio = 0.0;
        int numGibbsMoves = iterations;

        boolean[] tabu = new boolean[groupSizes.length];
        for (int i = 0; i < groupSizes.length; i++)
            tabu[i] = false;

        while (numGibbsMoves > 0) {

            int randomCoord = dr.math.MathUtils.nextInt(groupSizes.length);

            // do not bother to re-gibbs-sample coordinates when their
            // conditional distribution has not changed, i.e. when their
            // neighbours haven't been sampled
            if (!tabu[randomCoord]) {

                hastingsRatio += getSample(exponents, gammaRates, populationSizes, randomCoord);

                // make current index tabu, but un-tabuize its neighbours
                tabu[randomCoord] = true;
                if (randomCoord != 0)
                    tabu[randomCoord - 1] = false;
                if (randomCoord != groupSizes.length - 1)
                    tabu[randomCoord + 1] = false;
            }

            numGibbsMoves -= 1;

        }

        // store the new array of population sizes
        for (int j = 0; j < groupIndex; j++) {
            populationSizeParameter.setParameterValue(j, populationSizes[j]);
        }

        return hastingsRatio;

    }

    public final String getOperatorName() {
        return BayesianSkylineGibbsOperatorParser.BAYESIAN_SKYLINE_GIBBS_OPERATOR;
    }

}
