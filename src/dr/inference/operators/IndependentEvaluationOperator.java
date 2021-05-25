/*
 * GibbsIndependentNormalDistributionOperator.java
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

package dr.inference.operators;

import dr.inference.markovchain.Acceptor;
import dr.inference.markovchain.MarkovChain;
import dr.inference.mcmc.MCMC;
import dr.inference.mcmc.MCMCCriterion;
import dr.inference.model.CompoundLikelihood;
import dr.inference.model.Likelihood;
import dr.xml.*;

import java.util.ArrayList;
import java.util.List;

/**
 * This operator takes a list of likelihoods that can be updated independently using the
 * provided nested operator. It will attempt a Metropolis-Hastings update on each of the
 * likelihoods. It p a GibbsOperator as the calling MarkovChain assumes the acceptance
 * of the result.
 *
 * @author Andrew Rambaut
 *
 */
public class IndependentEvaluationOperator extends SimpleMCMCOperator implements GibbsOperator {

	public static final String OPERATOR_NAME = "independentEvaluationOperator";
	public static final String OPERATORS = "operators";

	private final CompoundLikelihood likelihoods;
	private final List<MCMCOperator> operators;
	private final Acceptor acceptor = new MCMCCriterion();

	public IndependentEvaluationOperator(CompoundLikelihood likelihoods, List<MCMCOperator> operators, double weight) {
		super(weight);

		this.likelihoods = likelihoods;
		this.operators = operators;
	}

	public String getPerformanceSuggestion() {
		return "";
	}

	public String getOperatorName() {
		return "IndependentEvaluationOperator(" + likelihoods.getId() + ")";
	}

	public int getStepCount() {
		return likelihoods.getLikelihoodCount();
	}

	/**
	 * change the parameter and return the hastings ratio.
	 */
	public double doOperation() {

		double[] oldLogLs = new double[likelihoods.getLikelihoodCount()];
		double[] hastingsRatios = new double[likelihoods.getLikelihoodCount()];

		// Get the likelihoods before the operator updates
		int i = 0;
		for (Likelihood likelihood: likelihoods.getLikelihoods()) {
			oldLogLs[i] = likelihood.getLogLikelihood();
			i++;
		}

		// Do an operation for each of the operators
		i = 0;
		for (MCMCOperator operator: operators) {
			hastingsRatios[i] = operator.operate();
			i++;
		}

		// Evaluate the updated states (in parallel - because that is the point of this
		likelihoods.getLogLikelihood();

		double[] logr = {0.0};

		// Decide whether to accept or reject each individual move
		i = 0;
		for (Likelihood likelihood: likelihoods.getLikelihoods()) {
			MCMCOperator operator = operators.get(i);

			double logL = likelihood.getLogLikelihood();
			boolean accept = acceptor.accept(oldLogLs[i], logL, hastingsRatios[i], logr);

			if (accept) {
				operator.accept(logL - oldLogLs[i]);
			} else {
				operator.reject();

				likelihood.getModel().restoreModelState();
			}

			if (operator instanceof AdaptableMCMCOperator) {
				adaptAcceptanceProbability((AdaptableMCMCOperator) operator, logr[0]);
			}

			i++;
		}

		return 0;
	}

	public void adaptAcceptanceProbability(AdaptableMCMCOperator op, double logr) {

		final boolean isAdaptable = op.getMode() == AdaptationMode.ADAPTATION_ON || (op.getMode() != AdaptationMode.ADAPTATION_OFF);

		if (isAdaptable) {
			final double p = op.getAdaptableParameter();

			final double i = OperatorSchedule.DEFAULT_TRANSFORM.transform(op.getAdaptationCount() + 2);

			double acceptance = Math.exp(logr);

			final double target = op.getTargetAcceptanceProbability();

			final double newp = p + ((1.0 / i) * (acceptance - target));

			if (newp > -Double.MAX_VALUE && newp < Double.MAX_VALUE) {
				op.setAdaptableParameter(newp);
			}
		}
	}

	private boolean isAdaptable(AdaptableMCMCOperator op) {
		return op.getMode() == AdaptationMode.ADAPTATION_ON
				|| (op.getMode() != AdaptationMode.ADAPTATION_OFF);
	}


	public static dr.xml.XMLObjectParser PARSER = new dr.xml.AbstractXMLObjectParser() {

		public String getParserName() {
			return OPERATOR_NAME;
		}

		public Object parseXMLObject(XMLObject xo) throws XMLParseException {

			double weight = xo.getDoubleAttribute(WEIGHT);

			CompoundLikelihood likelihoods = (CompoundLikelihood)xo.getChild(CompoundLikelihood.class);
			List<MCMCOperator> operators = new ArrayList<MCMCOperator>();
			XMLObject cxo = xo.getChild(OPERATORS);
			for (Object operator : cxo.getAllChildren(MCMCOperator.class)) {
				operators.add((MCMCOperator)operator);
			}

			if (likelihoods.getLikelihoodCount() != operators.size()) {
				throw new XMLParseException("The number of likelihoods does not match the number of operators");
			}

			return new IndependentEvaluationOperator(likelihoods, operators, weight);
		}

		//************************************************************************
		// AbstractXMLObjectParser implementation
		//************************************************************************

		public XMLSyntaxRule[] getSyntaxRules() {
			return rules;
		}

		private final XMLSyntaxRule[] rules = {
				AttributeRule.newDoubleRule(WEIGHT),
				new ElementRule(CompoundLikelihood.class),
				new ElementRule(OPERATORS, new XMLSyntaxRule[]{
						new ElementRule(MCMCOperator.class, 1, Integer.MAX_VALUE)
				}, "Operators to control"),
		};

		public String getParserDescription() {
			return "This element returns a sampler, disguised as a Gibbs operator, that controls operators on independent likelihoods.";
		}

		public Class getReturnType() {
			return MCMCOperator.class;
		}

	};

}
