package dr.evomodel.coalescent.operators;

import dr.inference.model.Parameter;
import dr.inference.operators.*;
import dr.evomodel.tree.TreeModel;
import dr.evomodel.coalescent.BayesianSkylineLikelihood;
import dr.evolution.coalescent.ConstantPopulation;
import dr.evolution.util.Units;
import dr.math.GammaDistribution;
import dr.xml.*;

/**
 * A Gibbs operator for the Bayesian Skyline model, inspired by the Hey &
 * Nielsen PNAS paper.
 * 
 * The operator is implemented as a standard proposal generator, returning a Hastings ratio.
 * The advantage of this scheme is that additional priors can be included without modifying 
 * the operator.  This is only be efficient for 'mild' priors, since acceptance rates 
 * easily get very low because this operator changes many parameters simultaneously.
 * 
 * The proposal generator can be modified by providing upper and lower bounds.
 * Four types of prior are implemented, so that the proposal behaves as a true Gibbs operator
 * if the chosen prior matches the actual prior.
 * 
 * jeffreys=	exponentialMarkov=
 * false	false			Uniform prior
 * true		false			Jeffrey's prior on the geometric average population size
 * false	true			An exponential Markov prior (as ExponentialMarkov)
 * true		true			Same, with a Jeffrey's prior on the first population size
 *  
 * @author Gerton Lunter
 */

public class BayesianSkylineGibbsOperator extends SimpleMCMCOperator {

	public static final String BAYESIAN_SKYLINE_GIBBS_OPERATOR = "generalizedSkylineGibbsOperator";

	public static final String POPULATION_SIZES = "populationSizes";

	public static final String GROUP_SIZES = "groupSizes";

	public static final String LOWER = "lower";

	public static final String UPPER = "upper";

	public static final String JEFFREYS = "Jeffreys";
	
	public static final String EXPONENTIALMARKOV = "exponentialMarkov";

	public static final String TYPE = "type";

	public static final String STEPWISE = "stepwise";

	public static final String LINEAR = "linear";

	public static final String EXPONENTIAL = "exponential";

	public static final int STEPWISE_TYPE = 0;

	public static final int LINEAR_TYPE = 1;

	public static final int EXPONENTIAL_TYPE = 2;

	private BayesianSkylineLikelihood bayesianSkylineLikelihood;

	private Parameter populationSizeParameter;

	private Parameter groupSizeParameter;

	private double upperBound;

	private double lowerBound;

	private boolean jeffreysGeometricPrior;        /* Jeffrey's prior on the geometric average */
	
	private boolean jeffreysPrior;                 /* Jeffrey's prior on the first Ne */
	
	private boolean exponentialMarkovPrior;        /* Exponential Markov prior */

	private double oldLogL, newLogL;
	
	public BayesianSkylineGibbsOperator(
			BayesianSkylineLikelihood bayesianSkylineLikelihood,
			Parameter populationSizeParameter,
			Parameter groupSizeParameter, int type, int weight,
			double lowerBound, double upperBound, boolean jeffreysPrior, boolean exponentialMarkov) {

		this.bayesianSkylineLikelihood = bayesianSkylineLikelihood;
		this.populationSizeParameter = populationSizeParameter;
		this.groupSizeParameter = groupSizeParameter;
		this.weight = weight;
		this.lowerBound = lowerBound;
		this.upperBound = upperBound;
		this.exponentialMarkovPrior = exponentialMarkov;
		if (exponentialMarkov) {
		    this.jeffreysPrior = jeffreysPrior;
		    this.jeffreysGeometricPrior = false;
		} else {
		    this.jeffreysGeometricPrior = jeffreysPrior;
		    this.jeffreysPrior = false;
		}

		assert( populationSizeParameter != null);
		assert( groupSizeParameter != null);
		assert( populationSizeParameter.getDimension() == groupSizeParameter.getDimension() );
		
		if (type != STEPWISE_TYPE) {
			throw new IllegalArgumentException("Should only get stepwise type here - sorry.");
		}
		
		System.out.println("Using a Bayesian skyline Gibbs operator (lo="+lowerBound+", hi="+upperBound+
				   ", Jeffreys="+jeffreysPrior+" exponentialMarkov="+exponentialMarkov+")");
	}

	
	public final String getPerformanceSuggestion() {
		return "This operator cannot be optimized";
	}

	
	double getSample( double gammaShape, double gammaRate, int groupIndex, double oldprevpopsize, double newprevpopsize, double popsize ) {
	    
	    // Add Jeffrey's prior, if required
	    if (jeffreysGeometricPrior) {
		gammaShape += 1.0 / groupSizeParameter.getDimension();
	    }
	    if (jeffreysPrior && groupIndex == 0) {
		gammaShape += 1.0;
	    }
			    
	    if (gammaShape < 0.0) {
		throw new Error("Bad shape parameter!  Bug in program!");
	    }

	    // Now sample from the Gamma distribution, or from the Gamma distribution with exponential markov prior
	    double sample, priordensity;
	    int iters = 200;
	    do {

		if (groupIndex > 0 && exponentialMarkovPrior) {
		    if (gammaShape < 1.0) {
			throw new IllegalArgumentException("Group size must be >= 2 when the exponential Markov prior is used.");
		    }
		    sample = 1.0 / GammaDistribution.nextExpGamma( gammaShape, 1.0 / gammaRate, newprevpopsize );
		} else {
		    sample = 1.0 / GammaDistribution.nextGamma( gammaShape, 1.0 / gammaRate );
		}
		iters -= 1;
		
	    } while ((sample < lowerBound || sample > upperBound) && iters>0);

	    if (iters==0) {
		// no sample within bounds found - fail
		return -1;
	    }
	    
	    // Calculate and update old and new log likelihoods
	    oldLogL += (-(gammaShape+1) * Math.log(popsize)) - gammaRate/popsize;
	    newLogL += (-(gammaShape+1) * Math.log(sample)) - gammaRate/sample;
	    if (groupIndex > 0 && exponentialMarkovPrior) {
		oldLogL += dr.math.ExponentialDistribution.logPdf(popsize, 1.0/oldprevpopsize);
		newLogL += dr.math.ExponentialDistribution.logPdf(sample, 1.0/newprevpopsize);
	    }
	    
	    //System.out.println(" i="+groupIndex+" shape="+gammaShape+" rate="+gammaRate+" old="+populationSizeParameter.getParameterValue(groupIndex)+" new="+sample);
	    return sample;
	}

	public double doOperation() throws OperatorFailedException {

		if (!bayesianSkylineLikelihood.getIntervalsKnown())
			bayesianSkylineLikelihood.setupIntervals();
				
		assert( populationSizeParameter != null);
		assert( groupSizeParameter != null);
		
		// Now enter a loop similar to the likelihood calculation, except that
		// we now just collect
		// waiting times, and exponents of the population size parameters, for
		// each group.

		oldLogL = 0.0;   // calculate the likelihood of the old and new state
		newLogL = 0.0;

		double oldprevpopsize = 0.0;      // passed to getSample to implement the exponential Markov prior
		double newprevpopsize = 0.0;
		double currentTime = 0.0;

		int groupIndex = 0;
		int[] groupSizes = bayesianSkylineLikelihood.getGroupSizes();
		double[] groupEnds = bayesianSkylineLikelihood.getGroupHeights();
		int subIndex = 0;

		ConstantPopulation cp = new ConstantPopulation(Units.Type.YEARS);
		double gammaShape = -1.0;    /* shape parameter -1.0 corresponds to exponent -2.0 in Gamma PDF, i.e. uniform on Ne */
		double gammaRate = 0.0;
		
		double[] newPopulationSizes = new double[bayesianSkylineLikelihood.getIntervalCount()];

		for (int j = 0; j < bayesianSkylineLikelihood.getIntervalCount(); j++) {

			// set the population size to the size of the middle of the current
			// interval
		    	double curpopsize = bayesianSkylineLikelihood.getPopSize(groupIndex,
		    		        	                              currentTime + (bayesianSkylineLikelihood.getInterval(j) / 2.0),
		    		        	                              groupEnds); 
			cp.setN0(curpopsize);

			gammaShape += bayesianSkylineLikelihood.calculateIntervalShapeParameter(cp,
							bayesianSkylineLikelihood.getInterval(j),
							currentTime, 
							bayesianSkylineLikelihood.getLineageCount(j),
							bayesianSkylineLikelihood.getIntervalType(j));
			gammaRate += bayesianSkylineLikelihood.calculateIntervalRateParameter(cp,
							bayesianSkylineLikelihood.getInterval(j),
							currentTime, 
							bayesianSkylineLikelihood.getLineageCount(j),
							bayesianSkylineLikelihood.getIntervalType(j)) * curpopsize;

			if (bayesianSkylineLikelihood.getIntervalType(j) == BayesianSkylineLikelihood.COALESCENT) {
			    subIndex += 1;
			    if (subIndex >= groupSizes[groupIndex]) {
				    
				// Finished with this group.  Compute a Gibbs sample
				double sample = getSample( gammaShape, gammaRate, groupIndex, oldprevpopsize, newprevpopsize, curpopsize );
				if (sample<0) throw new OperatorFailedException("Rejection sampling took too long.");

				newPopulationSizes[groupIndex] = sample;

				// Reset accumulators
				gammaShape = -1.0;
				gammaRate = 0.0;

				// Next group
				groupIndex += 1;
				subIndex = 0;
				oldprevpopsize = curpopsize;
				newprevpopsize = sample;
				
			    }
			}

			// insert zero-length coalescent intervals (for multifurcating
			// trees)
			int diff = bayesianSkylineLikelihood.getCoalescentEvents(j) - 1;
			for (int k = 0; k < diff; k++) {

			    	curpopsize = bayesianSkylineLikelihood.getPopSize(groupIndex,
			    						       currentTime + (bayesianSkylineLikelihood.getInterval(j) / 2.0),
			    						       groupEnds); 
				cp.setN0(curpopsize);

				gammaShape += bayesianSkylineLikelihood.calculateIntervalShapeParameter(cp,
					bayesianSkylineLikelihood.getInterval(j),
					currentTime, 
					bayesianSkylineLikelihood.getLineageCount(j),
					bayesianSkylineLikelihood.getIntervalType(j));
				gammaRate += bayesianSkylineLikelihood.calculateIntervalRateParameter(cp,
					bayesianSkylineLikelihood.getInterval(j),
					currentTime, 
					bayesianSkylineLikelihood.getLineageCount(j),
					bayesianSkylineLikelihood.getIntervalType(j)) * curpopsize;
			
				subIndex += 1;
				if (subIndex >= groupSizes[groupIndex]) {

					// Finished with this group.  Compute a Gibbs sample
					double sample = getSample( gammaShape, gammaRate, groupIndex, oldprevpopsize, newprevpopsize, curpopsize );
					if (sample<0) throw new OperatorFailedException("Rejection sampling took too long.");

					newPopulationSizes[groupIndex] = sample;

					// Reset accumulators
					gammaShape = -1.0;
					gammaRate = 0.0;

					// Next group
					groupIndex += 1;
					subIndex = 0;
					oldprevpopsize = curpopsize;
					newprevpopsize = sample;
				}

			}

			currentTime += bayesianSkylineLikelihood.getInterval(j);

		}

		for (int j = 0; j < groupIndex; j++) {
			populationSizeParameter.setParameterValue(j, newPopulationSizes[j]);
		}

		//System.out.println("Old + new likelihood = "+oldLogL+" "+newLogL);

		// Return a Hastings ratio equal to the ratio of old and new likelihoods, so that this move
		// will always be accepted (when priors match)
		return +(oldLogL - newLogL);
	}

	public final String getOperatorName() {
		return BAYESIAN_SKYLINE_GIBBS_OPERATOR;
	}

	public static dr.xml.XMLObjectParser PARSER = new dr.xml.AbstractXMLObjectParser() {

		public String getParserName() {
			return BAYESIAN_SKYLINE_GIBBS_OPERATOR;
		}

		public Object parseXMLObject(XMLObject xo) throws XMLParseException {

			double lowerBound = 0.0;
			double upperBound = Double.MAX_VALUE;
			boolean jeffreysPrior = true;
			boolean exponentialMarkovPrior = false;

			int weight = xo.getIntegerAttribute(WEIGHT);
			if (xo.hasAttribute(LOWER)) {
				lowerBound = xo.getDoubleAttribute(LOWER);
			}
			if (xo.hasAttribute(UPPER)) {
				upperBound = xo.getDoubleAttribute(UPPER);
			}
			if (xo.hasAttribute(JEFFREYS)) {
				jeffreysPrior = xo.getBooleanAttribute(JEFFREYS);
			}
			if (xo.hasAttribute(EXPONENTIALMARKOV)) {
				exponentialMarkovPrior = xo.getBooleanAttribute(EXPONENTIALMARKOV);
			}

			BayesianSkylineLikelihood bayesianSkylineLikelihood =(BayesianSkylineLikelihood) xo.getChild(BayesianSkylineLikelihood.class);

			// This is the parameter on which this operator acts
			Parameter paramPops = (Parameter) xo.getChild(Parameter.class);
			
			Parameter paramGroups = bayesianSkylineLikelihood.getGroupSizeParameter();

			int type = bayesianSkylineLikelihood.getType();

			if (type != BayesianSkylineLikelihood.STEPWISE_TYPE) {
				throw new XMLParseException("Need stepwise control points (set 'linear=\"false\"' in skyline Gibbs operator)");
			}

			return new BayesianSkylineGibbsOperator(bayesianSkylineLikelihood,
					paramPops, paramGroups, type, weight,
					lowerBound, upperBound, jeffreysPrior, exponentialMarkovPrior);

		}

		// ************************************************************************
		// AbstractXMLObjectParser implementation
		// ************************************************************************

		public String getParserDescription() {
			return "This element returns a Gibbs operator for the joint distribution of population sizes.";
		}

		public Class getReturnType() {
			return MCMCOperator.class;
		}

		public XMLSyntaxRule[] getSyntaxRules() {
			return rules;
		}

		private XMLSyntaxRule[] rules = new XMLSyntaxRule[] {
				AttributeRule.newBooleanRule(LINEAR, true),
				AttributeRule.newIntegerRule(WEIGHT),
				AttributeRule.newDoubleRule(LOWER),
				AttributeRule.newDoubleRule(UPPER),
				AttributeRule.newBooleanRule(JEFFREYS, true),
				AttributeRule.newBooleanRule(EXPONENTIALMARKOV,true),
				new ElementRule(BayesianSkylineLikelihood.class),
				new ElementRule(Parameter.class)
			};

	};

}
