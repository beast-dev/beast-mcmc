/*
 * ParallelCompoundLikelihood.java
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

package dr.inference.model;

import dr.inference.parallel.MPIServices;
import dr.xml.*;

import java.util.*;

/**
 * @author Marc A. Suchard
 */

public class ParallelCompoundLikelihood extends CompoundLikelihood {

	public static final String PARALLEL_COMPOUND_LIKELIHOOD = "parallelCompoundLikelihood";
	public static final String LOCAL_CHECK = "doLocalCheck";
	public static final String RUN_PARALLEL = "doInParallel";

	public ParallelCompoundLikelihood(Collection<Likelihood> likelihoods, boolean doParallel, boolean checkLocal) {
		super(1, likelihoods);
		this.doParallel = doParallel;
		this.checkLocal = checkLocal;
	}

	private boolean doParallel = true;
	private boolean checkLocal = false;

	public double getLogLikelihood() {

		double logLikelihood = 0;
		if (doParallel) {

			logLikelihood = getLogLikelihoodRemote();
			if (checkLocal) {
				super.makeDirty();
				double logLikelihoodLocal = super.getLogLikelihood();
				System.err.printf("Local: %5.4f  Remote: %5.4f\n", logLikelihoodLocal, logLikelihood);
			}

		} else
			logLikelihood = super.getLogLikelihood();

		return logLikelihood;
	}

	private double getLogLikelihoodRemote() {
		double logLikelihood = 0.0;

		final int N = getLikelihoodCount();

		List<ParallelLikelihood> likelihoodsDistributed = new ArrayList<ParallelLikelihood>();
		List<Integer> processorList = new ArrayList<Integer>();

		for (int i = 0; i < N; i++) {
			ParallelLikelihood likelihood = (ParallelLikelihood) getLikelihood(i);
			if (!likelihood.getLikelihoodKnown()) {
				//    if (true) {
				final int processor = i + 1;
//				MPIServices.requestLikelihood(processor);
//				((AbstractModel) getLikelihood(i).getModel()).sendState(processor);
				likelihoodsDistributed.add(likelihood);
				processorList.add(processor);
			} else {
				logLikelihood += likelihood.getLogLikelihood();
			}
		}

		final int size = likelihoodsDistributed.size();

		if (size == 1) { // only one, so do locally

			logLikelihood += likelihoodsDistributed.get(0).getLogLikelihood();

		} else if (size > 1) {

			// Distribute calculations
			int index = 0;
			for (ParallelLikelihood likelihood : likelihoodsDistributed) {
				int processor = processorList.get(index++);
				MPIServices.requestLikelihood(processor);
				((AbstractModel) likelihood.getModel()).sendState(processor);
			}

			// Implicit barrier

			// Collect calculations
			index = 0;
			for (ParallelLikelihood likelihood : likelihoodsDistributed) {
				int processor = processorList.get(index++);
				double l = MPIServices.receiveDouble(processor);
				logLikelihood += l;
				likelihood.setLikelihood(l);         // todo don't we need to set all of the submodels ????
			}

			// todo Use Gather instead of blocking loop

		}

		return logLikelihood;
	}


	public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

		public String getParserName() {
			return PARALLEL_COMPOUND_LIKELIHOOD;
		}

		//      public String[] getParserNames() { return new String[] { getParserName(), "posterior", "prior", "likelihood" }; }

		public Object parseXMLObject(XMLObject xo) throws XMLParseException {

			boolean doParallel = true;
			boolean checkLocal = false;

			if (xo.hasAttribute(LOCAL_CHECK)) {
				checkLocal = xo.getBooleanAttribute(LOCAL_CHECK);
			}

			if (xo.hasAttribute(RUN_PARALLEL)) {
				doParallel = xo.getBooleanAttribute(RUN_PARALLEL);
			}
            List<Likelihood> likelihoods = new ArrayList<Likelihood>();
            for (int i = 0; i < xo.getChildCount(); i++) {
                if (xo.getChild(i) instanceof Likelihood) {
                    likelihoods.add((Likelihood) xo.getChild(i));
                } else {
                    throw new XMLParseException("An element which is not a likelihood has been added to a " + PARALLEL_COMPOUND_LIKELIHOOD + " element");
                }
            }

			ParallelCompoundLikelihood compoundLikelihood = new ParallelCompoundLikelihood(likelihoods, doParallel, checkLocal);

			return compoundLikelihood;
		}

		//************************************************************************
		// AbstractXMLObjectParser implementation
		//************************************************************************

		public String getParserDescription() {
			return "A likelihood function which is simply the product of its component likelihood functions.";
		}

		public XMLSyntaxRule[] getSyntaxRules() {
			return rules;
		}

		private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
				new ElementRule(Likelihood.class, 1, Integer.MAX_VALUE),
		};

		public Class getReturnType() {
			return ParallelCompoundLikelihood.class;
		}
	};
}
