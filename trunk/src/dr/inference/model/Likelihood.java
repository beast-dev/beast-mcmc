/*
 * Likelihood.java
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

import dr.inference.loggers.Loggable;
import dr.util.Identifiable;

/**
 * classes that calculate likelihoods should implement this interface.
 *
 * @author Alexei Drummond
 * @author Andrew Rambaut
 *
 * @version $Id: Likelihood.java,v 1.16 2005/05/24 20:26:00 rambaut Exp $
 */

public interface Likelihood extends Loggable, Identifiable {

	/**
	 * Get the model.
	 * @return the model.
	 */
	Model getModel();

	/**
	 * Get the log likelihood.
	 * @return the log likelihood.
	 */
	double getLogLikelihood();

	/**
	 * Forces a complete recalculation of the likelihood next time getLikelihood is called
	 */
	void makeDirty();

	/**
	 * A simple abstract base class for likelihood functions
	 */

	public abstract class Abstract implements Likelihood, ModelListener {

		public Abstract(Model model) {

			this.model = model;
			if (model != null) model.addModelListener(this);
		}

		public void modelChangedEvent(Model model, Object object, int index) {
			makeDirty();
		}

		// **************************************************************
	    // Likelihood IMPLEMENTATION
	    // **************************************************************

		/**
		 * Get the model.
		 * @return the model.
		 */
		public Model getModel() { return model; }

		public final double getLogLikelihood() {
			if (!getLikelihoodKnown()) {
				logLikelihood = calculateLogLikelihood();
				likelihoodKnown = true;
			}
			return logLikelihood;
		}

		public void makeDirty() {
			likelihoodKnown = false;
		}

		/**
		 * Called to decide if the likelihood must be calculated. Can be overridden
		 * (for example, to always return false).
         * @return  true if no need to recompute likelihood
         */
		protected boolean getLikelihoodKnown() {
			return likelihoodKnown;
		}

		protected abstract double calculateLogLikelihood();

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

		private Model model;
		private double logLikelihood;
		private boolean likelihoodKnown = false;
	}
}
