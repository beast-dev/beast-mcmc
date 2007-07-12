/*
 * CompoundLikelihood.java
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

package dr.inference.model;

import dr.xml.*;
import dr.util.NumberFormatter;

import java.util.ArrayList;
import java.util.List;
import java.util.logging.Logger;
import java.util.concurrent.*;

/**
 * A likelihood function which is simply the product of a set of likelihood functions.
 *
 * @author Alexei Drummond
 * @author Andrew Rambaut
 *
 * @version $Id: CompoundLikelihood.java,v 1.19 2005/05/25 09:14:36 rambaut Exp $
 */
public class CompoundLikelihood implements Likelihood {

	public static final String COMPOUND_LIKELIHOOD = "compoundLikelihood";
	public static final String THREADS = "threads";
	public static final String POSTERIOR = "posterior";
	public static final String PRIOR = "prior";
	public static final String LIKELIHOOD = "likelihood";

	public CompoundLikelihood(int threads) {
		if (threads > 1) {
			pool = Executors.newFixedThreadPool(threads);
//			pool = Executors.newCachedThreadPool();
		} else {
			pool = null;
		}
	}

	public void addLikelihood(Likelihood likelihood) {

		if ( !likelihoods.contains(likelihood) ) {

			likelihoods.add(likelihood);
			if (likelihood.getModel() != null) {
				compoundModel.addModel(likelihood.getModel());
			}

			likelihoodCallers.add(new LikelihoodCaller(likelihood));
		}
	}

	public int getLikelihoodCount() {
		return likelihoods.size();
	}

	public final Likelihood getLikelihood(int i) {
		return likelihoods.get(i);
	}

	// **************************************************************
	// Likelihood IMPLEMENTATION
	// **************************************************************

	public Model getModel() { return compoundModel; }

	public final double getLogLikelihood() {
		double logLikelihood = 0.0;

        if (pool == null) {
			// Single threaded

			//System.err.println("mixed of " + likelihoods.size());
			for (Likelihood likelihood : likelihoods) {
				double l = likelihood.getLogLikelihood();

				// if the likelihood is zero then short cut the rest of the likelihoods
				// This means that expensive likelihoods such as TreeLikelihoods should
				// be put after cheap ones such as BooleanLikelihoods
				if (l == Double.NEGATIVE_INFINITY) return Double.NEGATIVE_INFINITY;
				logLikelihood += l;
			}
		} else {
			try {
				List<Future<Double>> results = pool.invokeAll(likelihoodCallers);

                for (Future<Double> result : results) {

                    double logL = result.get();
                    logLikelihood += logL;
                }

			} catch (InterruptedException e) {
				e.printStackTrace();
			} catch (ExecutionException e) {
				e.printStackTrace();
			}
		}

        return logLikelihood;
	}

	public void makeDirty() {

		for (Likelihood likelihood : likelihoods) {
			likelihood.makeDirty();
		}
	}

	public String getDiagnosis() {
		String message = "";
		boolean first = true;

		for (Likelihood lik : likelihoods) {

			if (!first) {
				message += ", ";
			} else {
				first = false;
			}

			String id = lik.getId();
			if (id == null || id.trim().length() == 0) {
				String[] parts = lik.getClass().getName().split("\\.");
				id = parts[parts.length - 1];
			}

			message += id + "=";


			if (lik instanceof CompoundLikelihood) {
				String d = ((CompoundLikelihood)lik).getDiagnosis();
				if (d != null && d.length() > 0) {
					message += "(" + d + ")";
				}
			} else {

				if (lik.getLogLikelihood() == Double.NEGATIVE_INFINITY) {
					message += "-Inf";
				} else if (Double.isNaN(lik.getLogLikelihood())) {
					message += "NaN";
				} else {
					NumberFormatter nf = new NumberFormatter(6);
					message += nf.formatDecimal(lik.getLogLikelihood(), 4);
				}
			}
		}

		return message;
	}

	public String toString() {

		return Double.toString(getLogLikelihood());

	}

	// **************************************************************
	// Loggable IMPLEMENTATION
	// **************************************************************

	/**
	 * @return the log columns.
	 */
	public dr.inference.loggers.LogColumn[] getColumns() {
		return new dr.inference.loggers.LogColumn[] {
				new LikelihoodColumn(getId())
		};
	}

	private class LikelihoodColumn extends dr.inference.loggers.NumberColumn {
		public LikelihoodColumn(String label) { super(label); }
		public double getDoubleValue() { return getLogLikelihood(); }
	}

	// **************************************************************
	// Identifiable IMPLEMENTATION
	// **************************************************************

	private String id = null;

	public void setId(String id) { this.id = id; }

	public String getId() { return id; }

	public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

		public String getParserName() { return COMPOUND_LIKELIHOOD; }

		public String[] getParserNames() { return new String[] { getParserName(), "posterior", "prior", "likelihood" }; }

		public Object parseXMLObject(XMLObject xo) throws XMLParseException {

			int threads = 0;

			if (xo.hasAttribute(THREADS)) {
				threads = xo.getIntegerAttribute(THREADS);
			}
			CompoundLikelihood compoundLikelihood = new CompoundLikelihood(threads);

			for (int i = 0; i < xo.getChildCount(); i++) {
				if (xo.getChild(i) instanceof Likelihood) {
					compoundLikelihood.addLikelihood((Likelihood)xo.getChild(i));
				} else {

					Object rogueElement = xo.getChild(i);

					throw new XMLParseException("An element (" + rogueElement + ") which is not a likelihood has been added to a " + COMPOUND_LIKELIHOOD + " element");
				}
			}

            if (threads > 1) {
                Logger.getLogger("dr.evomodel").info("Likelihood is using " + threads + " threads.");
            }


            return compoundLikelihood;
		}

		//************************************************************************
		// AbstractXMLObjectParser implementation
		//************************************************************************

		public String getParserDescription() {
			return "A likelihood function which is simply the product of its component likelihood functions.";
		}

		public XMLSyntaxRule[] getSyntaxRules() { return rules; }

		private XMLSyntaxRule[] rules = new XMLSyntaxRule[] {
				new ElementRule(Likelihood.class, 1, Integer.MAX_VALUE ),
				AttributeRule.newIntegerRule(THREADS, true)
		};

		public Class getReturnType() { return CompoundLikelihood.class; }
	};

	private final ExecutorService pool;

	private ArrayList<Likelihood> likelihoods = new ArrayList<Likelihood>();
	private CompoundModel compoundModel = new CompoundModel("compoundModel");

	private List<Callable<Double>> likelihoodCallers = new ArrayList<Callable<Double>>();
	class LikelihoodCaller implements Callable<Double> {

		public LikelihoodCaller(Likelihood likelihood) {
			this.likelihood = likelihood;
		}

		public Double call() throws Exception {
			return likelihood.getLogLikelihood();
		}

		private final Likelihood likelihood;
	}
}

