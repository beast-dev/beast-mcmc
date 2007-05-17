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
 * Note: If the uniform prior is used ('jeffreys="false"), then the minimum
 * group size is 2; sampling breaks down for groups containing just 1 coalescent
 * (although with bounds, this could be repaired).
 * 
 * Note: Currently, the tree must be provided, as well as the Bayesian Skyline
 * likelihood object. This might change.
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

	private boolean jeffreysPrior;

	public BayesianSkylineGibbsOperator(
			BayesianSkylineLikelihood bayesianSkylineLikelihood,
			Parameter populationSizeParameter,
			Parameter groupSizeParameter, int type, int weight,
			double lowerBound, double upperBound, boolean jeffreysPrior) {

		this.bayesianSkylineLikelihood = bayesianSkylineLikelihood;
		this.populationSizeParameter = populationSizeParameter;
		this.groupSizeParameter = groupSizeParameter;
		this.weight = weight;
		this.lowerBound = lowerBound;
		this.upperBound = upperBound;
		this.jeffreysPrior = jeffreysPrior;

		assert( populationSizeParameter != null);
		assert( groupSizeParameter != null);
		assert( populationSizeParameter.getDimension() == groupSizeParameter.getDimension() );
		
		if (type != STEPWISE_TYPE) {
			throw new IllegalArgumentException(
					"Should only get stepwise type here - sorry.");
		}
		
		System.out.println("Using a Bayesian skyline Gibbs operator (lo="+lowerBound+", hi="+upperBound+", Jeffrey's="+jeffreysPrior+")");

	}

	public final String getPerformanceSuggestion() {
		return "Always use Gibbs moves";
	}

	double getSample( double gammaShape, double gammaRate ) {
	    
	    // Add Jeffrey's prior, if required
	    if (jeffreysPrior) {
		gammaShape += 1.0 / groupSizeParameter.getDimension();
	    }
	    
	    if (gammaShape < 0.0) {
		throw new Error("For BayesianSkylineGibbs with uniform priors, minimum group size is 2");
	    }

	    // Now sample from the Gamma distribution
	    double sample;
	    int iters = 500;
	    do {
		sample = 1.0 / GammaDistribution.nextGamma( gammaShape, 1.0 / gammaRate );
		iters -= 1;
	    } while ((sample < lowerBound || sample > upperBound) && iters>0);

	    // This is a hack -- really, I should return a sample from the extreme value distribution.
	    if (sample < lowerBound) sample = lowerBound;
	    if (sample > upperBound) sample = upperBound;

	    //System.out.println(" i="+j+" shape="+gammaShape+" rate="+gammaRate+" old="+populationSizeParameter.getParameterValue(groupIndex)+" new="+sample);
	    return sample;
	}

	public double doOperation() throws OperatorFailedException {

		if (!bayesianSkylineLikelihood.getIntervalsKnown())
			bayesianSkylineLikelihood.setupIntervals();
		
		assert( populationSizeParameter != null);
		assert( groupSizeParameter != null);
		//System.out.println("Old likelihood = "+bayesianSkylineLikelihood.calculateLogLikelihood());

		// Now enter a loop similar to the likelihood calculation, except that
		// we now just collect
		// waiting times, and exponents of the population size parameters, for
		// each group.

		double oldLogL = 0.0;   // calculate the likelihood of the old and new state
		double newLogL = 0.0;
		
		double currentTime = 0.0;

		int groupIndex = 0;
		int[] groupSizes = bayesianSkylineLikelihood.getGroupSizes();
		double[] groupEnds = bayesianSkylineLikelihood.getGroupHeights();
		int subIndex = 0;

		ConstantPopulation cp = new ConstantPopulation(Units.Type.YEARS);
		double gammaShape = -1.0;    /* shape parameter -1.0 corresponds to exponent -2.0 in Gamma PDF, i.e. uniform on Ne */
		double gammaRate = 0.0;

		for (int j = 0; j < bayesianSkylineLikelihood.getIntervalCount(); j++) {

			// set the population size to the size of the middle of the current
			// interval
		    	double popsize = bayesianSkylineLikelihood.getPopSize(groupIndex,
		    		        	                              currentTime + (bayesianSkylineLikelihood.getInterval(j) / 2.0),
		    		        	                              groupEnds); 
			cp.setN0(popsize);

			gammaShape += bayesianSkylineLikelihood.calculateIntervalShapeParameter(cp,
							bayesianSkylineLikelihood.getInterval(j),
							currentTime, 
							bayesianSkylineLikelihood.getLineageCount(j),
							bayesianSkylineLikelihood.getIntervalType(j));
			gammaRate += bayesianSkylineLikelihood.calculateIntervalRateParameter(cp,
							bayesianSkylineLikelihood.getInterval(j),
							currentTime, 
							bayesianSkylineLikelihood.getLineageCount(j),
							bayesianSkylineLikelihood.getIntervalType(j)) * popsize;

			if (bayesianSkylineLikelihood.getIntervalType(j) == BayesianSkylineLikelihood.COALESCENT) {
			    subIndex += 1;
			    if (subIndex >= groupSizes[groupIndex]) {
				    
				// Finished with this group.  Compute a Gibbs sample
				double sample = getSample( gammaShape, gammaRate );

				// Calculate old and new log likelihoods
				oldLogL += (-(gammaShape+1) * Math.log(popsize)) - gammaRate/popsize;
				newLogL += (-(gammaShape+1) * Math.log(sample)) - gammaRate/sample;

				populationSizeParameter.setParameterValue(groupIndex, sample);

				// Reset accumulators
				gammaShape = -1.0;
				gammaRate = 0.0;

				// Next group
				groupIndex += 1;
				subIndex = 0;
				
			    }
			}

			// insert zero-length coalescent intervals (for multifurcating
			// trees)
			int diff = bayesianSkylineLikelihood.getCoalescentEvents(j) - 1;
			for (int k = 0; k < diff; k++) {

			    	popsize = bayesianSkylineLikelihood.getPopSize(groupIndex,
			    						       currentTime + (bayesianSkylineLikelihood.getInterval(j) / 2.0),
			    						       groupEnds); 
				cp.setN0(popsize);

				gammaShape += bayesianSkylineLikelihood.calculateIntervalShapeParameter(cp,
					bayesianSkylineLikelihood.getInterval(j),
					currentTime, 
					bayesianSkylineLikelihood.getLineageCount(j),
					bayesianSkylineLikelihood.getIntervalType(j));
				gammaRate += bayesianSkylineLikelihood.calculateIntervalRateParameter(cp,
					bayesianSkylineLikelihood.getInterval(j),
					currentTime, 
					bayesianSkylineLikelihood.getLineageCount(j),
					bayesianSkylineLikelihood.getIntervalType(j)) * popsize;
			
				subIndex += 1;
				if (subIndex >= groupSizes[groupIndex]) {

					// Finished with this group.  Compute a Gibbs sample
					double sample = getSample( gammaShape, gammaRate );

					// Calculate old and new log likelihoods
					oldLogL += (-(gammaShape+1) * Math.log(popsize)) - gammaRate/popsize;
					newLogL += (-(gammaShape+1) * Math.log(sample)) - gammaRate/sample;
					
					populationSizeParameter.setParameterValue(groupIndex, sample);

					// Reset accumulators
					gammaShape = -1.0;
					gammaRate = 0.0;

					// Next group
					groupIndex += 1;
					subIndex = 0;
				}

			}

			currentTime += bayesianSkylineLikelihood.getInterval(j);

		}

		//System.out.println("Old + new likelihood = "+oldLogL+" "+newLogL);

		// Return a Hastings ratio equal to the ratio of old and new likelihoods, so that this move
		// will always be accepted (if no priors on Ne have been defined)
		return oldLogL - newLogL;
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

			BayesianSkylineLikelihood bayesianSkylineLikelihood =(BayesianSkylineLikelihood) xo.getChild(BayesianSkylineLikelihood.class);

			// This is the parameter on which this operator acts
			Parameter paramPops = (Parameter) xo.getChild(Parameter.class);
			
			Parameter paramGroups = bayesianSkylineLikelihood.getGroupSizeParameter();

			int type = bayesianSkylineLikelihood.getType();

			if (type != BayesianSkylineLikelihood.STEPWISE_TYPE) {
				throw new XMLParseException(
						"Need stepwise control points (set 'linear=\"false\"' in skyline Gibbs operator)");
			}

			return new BayesianSkylineGibbsOperator(bayesianSkylineLikelihood,
					paramPops, paramGroups, type, weight,
					lowerBound, upperBound, jeffreysPrior);

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
				new ElementRule(BayesianSkylineLikelihood.class),
				new ElementRule(Parameter.class)
			};

	};

}
