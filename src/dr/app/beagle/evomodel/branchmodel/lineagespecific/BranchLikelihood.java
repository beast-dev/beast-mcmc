/*
 * BranchLikelihood.java
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

package dr.app.beagle.evomodel.branchmodel.lineagespecific;

import java.util.LinkedList;
import java.util.List;

import dr.evolution.tree.NodeRef;
import dr.evomodel.branchratemodel.CountableBranchCategoryProvider;
import dr.evomodel.tree.TreeModel;
import dr.inference.loggers.LogColumn;
import dr.inference.loggers.NumberColumn;
import dr.inference.model.CompoundModel;
import dr.inference.model.Likelihood;
import dr.inference.model.Model;
import dr.inference.model.Parameter;

@SuppressWarnings("serial")
public class BranchLikelihood implements Likelihood {

	private TreeModel treeModel;
	private List<Likelihood> uniqueLikelihoods;

	private List<Likelihood> branchLikelihoods;
	private Parameter categoriesParameter;

	private String id = null;
	private boolean used = false;

	// for discrete categories
	private CountableBranchCategoryProvider categoriesProvider;

	private final CompoundModel compoundModel = new CompoundModel(
			"compoundModel");

	public BranchLikelihood(TreeModel treeModel, List<Likelihood> likelihoods,
			Parameter categoriesParameter) {

		this.treeModel = treeModel;
		this.uniqueLikelihoods = likelihoods;
		this.categoriesParameter = categoriesParameter;

		if (this.treeModel != null) {

			this.categoriesProvider = new CountableBranchCategoryProvider.IndependentBranchCategoryModel(
					treeModel, categoriesParameter);

		}

		this.branchLikelihoods = getBranchLikelihoods();

	}// END: Constructor

	public final Likelihood getLikelihood(int index) {
		return branchLikelihoods.get(index);
	}
	
	private List<Likelihood> getBranchLikelihoods() {

		// linked list preserves order
		List<Likelihood> loglikes = new LinkedList<Likelihood>();

		if (treeModel != null) {

			for (NodeRef branch : treeModel.getNodes()) {

				if (!treeModel.isRoot(branch)) {

					int branchCategory = categoriesProvider.getBranchCategory(
							treeModel, branch);
					int index = (int) categoriesParameter
							.getParameterValue(branchCategory);

					Likelihood branchLikelihood = uniqueLikelihoods.get(index);

					loglikes.add(branchLikelihood);
					compoundModel.addModel(branchLikelihood.getModel());

				}
			}// END: branch loop

		} else {// if no tree then read them in supplied order

			int dim = categoriesParameter.getDimension();
			if (dim != uniqueLikelihoods.size()) {
				throw new RuntimeException("Dimensionality mismatch!");
			}// END: size of categoriesParameter check

			loglikes.addAll(uniqueLikelihoods);

		}// END: tree check

		return loglikes;
	}// END: getBranchLikelihoods

	@Override
	public double getLogLikelihood() {

		double loglike = 0;
		for (Likelihood like : branchLikelihoods) {
			loglike += like.getLogLikelihood();
		}

		return loglike;
	}// END: getLogLikelihood

	public int getLikelihoodCount() {
		return branchLikelihoods.size();
	}

	public List<Likelihood> getLikelihoods() {
		return branchLikelihoods;
	}

	@Override
	public LogColumn[] getColumns() {
		return new dr.inference.loggers.LogColumn[] { new LikelihoodColumn(
				getId() == null ? "likelihood" : getId()) };
	}// END: getColumns

	@Override
	public String getId() {
		return this.id;
	}// END: getId

	@Override
	public void setId(String id) {
		this.id = id;
	}// END: setId

	@Override
	public Model getModel() {
		return compoundModel;
	}

	@Override
	public void makeDirty() {

		for (Likelihood likelihood : uniqueLikelihoods) {
			likelihood.makeDirty();
		}

	}// END: makeDirty

	@Override
	public String prettyName() {
		return Abstract.getPrettyName(this);
	}// END: prettyName

	@Override
	public boolean isUsed() {
		return used;
	}// END: isUsed

	@Override
	public void setUsed() {
		used = true;
		for (Likelihood like : branchLikelihoods) {
			like.setUsed();
		}

	}// END: setUsed

	@Override
	public boolean evaluateEarly() {
		return false;
	}

	private class LikelihoodColumn extends NumberColumn {

		public LikelihoodColumn(String label) {
			super(label);
		}// END: Constructor

		public double getDoubleValue() {
			return getLogLikelihood();
		}

	}// END: LikelihoodColumn class

}// END: class
