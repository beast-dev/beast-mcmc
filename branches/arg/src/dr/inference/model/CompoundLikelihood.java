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

import dr.util.NumberFormatter;
import dr.xml.*;

import java.util.ArrayList;

/**
 * A likelihood function which is simply the product of a set of likelihood functions.
 *
 * @author Alexei Drummond
 * @author Andrew Rambaut
 * @version $Id: CompoundLikelihood.java,v 1.19 2005/05/25 09:14:36 rambaut Exp $
 */
public class CompoundLikelihood implements Likelihood {

	public static final String COMPOUND_LIKELIHOOD = "compoundLikelihood";
	public static final String POSTERIOR = "posterior";
	public static final String PRIOR = "prior";
	public static final String LIKELIHOOD = "likelihood";

	public CompoundLikelihood() {
	}

	public void addLikelihood(Likelihood likelihood) {

		if (!likelihoods.contains(likelihood)) {

			likelihoods.add(likelihood);
			if (likelihood.getModel() != null) {
				compoundModel.addModel(likelihood.getModel());
			}

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

	public Model getModel() {
		return compoundModel;
	}

	public double getLogLikelihood() {
		double logLikelihood = 0.0;

		for (int i = 0; i < likelihoods.size(); i++) {
			double l = getLikelihood(i).getLogLikelihood();

			// if the likelihood is zero then short cut the rest of the likelihoods
			// This means that expensive likelihoods such as TreeLikelihoods should
			// be put after cheap ones such as BooleanLikelihoods
			if (l == Double.NEGATIVE_INFINITY) return Double.NEGATIVE_INFINITY;

			logLikelihood += l;
		}

		return logLikelihood;
	}

	public void makeDirty() {

		for (int i = 0; i < likelihoods.size(); i++) {
			getLikelihood(i).makeDirty();
		}
	}

	public String getDiagnosis() {
		String message = "";
		boolean first = true;

		for (int i = 0; i < likelihoods.size(); i++) {
			Likelihood lik = getLikelihood(i);

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
				String d = ((CompoundLikelihood) lik).getDiagnosis();
				if (d != null & d.length() > 0) {
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
		return new dr.inference.loggers.LogColumn[]{
				new LikelihoodColumn(getId())
		};
	}

	private class LikelihoodColumn extends dr.inference.loggers.NumberColumn {
		public LikelihoodColumn(String label) {
			super(label);
		}

		public double getDoubleValue() {
			return getLogLikelihood();
		}
	}

	// **************************************************************
	// Identifiable IMPLEMENTATION
	// **************************************************************

	private String id = null;

	public void setId(String id) {
		this.id = id;
	}

	public String getId() {
		return id;
	}

	public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

		public String getParserName() {
			return COMPOUND_LIKELIHOOD;
		}

		public String[] getParserNames() {
			return new String[]{getParserName(), "posterior", "prior", "likelihood"};
		}

		public Object parseXMLObject(XMLObject xo) throws XMLParseException {

			CompoundLikelihood compoundLikelihood = new CompoundLikelihood();

			for (int i = 0; i < xo.getChildCount(); i++) {
				if (xo.getChild(i) instanceof Likelihood) {
					compoundLikelihood.addLikelihood((Likelihood) xo.getChild(i));
				} else {
					throw new XMLParseException("An element which is not a likelihood has been added to a " + COMPOUND_LIKELIHOOD + " element");
				}
			}
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
			return CompoundLikelihood.class;
		}
	};

	private ArrayList<Likelihood> likelihoods = new ArrayList<Likelihood>();
	private CompoundModel compoundModel = new CompoundModel("compoundModel");
}

