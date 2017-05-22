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

package dr.evomodel.branchmodel.lineagespecific;

import java.util.*;

import dr.evomodel.treelikelihood.BeagleTreeLikelihood;
import dr.evomodel.branchratemodel.CountableBranchCategoryProvider;
import dr.evomodel.tree.TreeModel;
import dr.inference.loggers.LogColumn;
import dr.inference.loggers.NumberColumn;
import dr.inference.model.CompoundModel;
import dr.inference.model.CompoundParameter;
import dr.inference.model.Likelihood;
import dr.inference.model.Model;
import dr.inference.model.Parameter;

@SuppressWarnings("serial")
public class BranchLikelihood implements Likelihood {

	// Constructor fields
//	private PatternList patternList;
	private TreeModel treeModel;
//	private BranchModel branchModel;
//	private SiteRateModel siteRateModel;
//	private FrequencyModel freqModel;
//	private BranchRateModel branchRateModel;
	
    private BeagleTreeLikelihood likelihood;
	private CountableBranchCategoryProvider categoriesProvider;

	private Parameter categoriesParameter;
	private CompoundParameter uniqueParameters;
	private CountableRealizationsParameter allParameters;
	
	private final CompoundModel compoundModel = new CompoundModel(
			"compoundModel");

	public BranchLikelihood(	
//			 PatternList patternList, //
			 TreeModel treeModel, //
//			 BranchModel branchModel, //
//			 SiteRateModel siteRateModel, //
//			 FrequencyModel freqModel, //
//			 BranchRateModel branchRateModel, //
			BeagleTreeLikelihood likelihood,
			 Parameter categoriesParameter,
			 CompoundParameter uniqueParameters,
			 CountableRealizationsParameter allParameters
			 ) {

//		this.patternList = patternList;
//		this.treeModel = treeModel;
//		this.branchModel = branchModel;
//		this.siteRateModel = siteRateModel;
//		this.freqModel = freqModel;
//		this.branchRateModel = branchRateModel;
		this.treeModel = treeModel;
		this.categoriesParameter = categoriesParameter;
		this.uniqueParameters = uniqueParameters;
		this.allParameters = allParameters;
		
			this.categoriesProvider = new CountableBranchCategoryProvider.IndependentBranchCategoryModel(
					this.treeModel, categoriesParameter);

	}// END: Constructor

	
	

	@Override
	public double getLogLikelihood() {
		
		double loglike = 0;

		for (int i = 0; i < categoriesParameter.getDimension(); i++) {
			
			int category = (int) categoriesParameter.getParameterValue(i);			
			double value = uniqueParameters.getParameterValue(category);
			
             allParameters.setParameterValue(i, value);
			
		}//END: i loop

		likelihood.makeDirty();
        loglike = likelihood.getLogLikelihood();

		return loglike;
	}// END: getLogLikelihood

	public double getLogLikelihood(int index) {
		
		
//		likelihood.getBranchModel().getSubstitutionModels().get(index).getVariable(0).setValue(index, value)
		
		double tmp1 = likelihood.getLogLikelihood();
		
		int category = (int) categoriesParameter.getParameterValue(index);			
		double value = uniqueParameters.getParameterValue(category);
		allParameters.setParameterValue(index, value);
		
		likelihood.makeDirty();
		
		
		double tmp2 = likelihood.getLogLikelihood();
		
		
		return tmp2-tmp1;
	}

	@Override
	public LogColumn[] getColumns() {
		return new dr.inference.loggers.LogColumn[] { new LikelihoodColumn(
				getId() == null ? "likelihood" : getId()) };
	}// END: getColumns

	@Override
	public String getId() {
		return"";
	}// END: getId

	@Override
	public void setId(String id) {
	}// END: setId

	@Override
	public Model getModel() {
		return compoundModel;
	}

	@Override
	public void makeDirty() {


	}// END: makeDirty

	@Override
	public String prettyName() {
		return Abstract.getPrettyName(this);
	}// END: prettyName

	@Override
	public Set<Likelihood> getLikelihoodSet() {
		return new HashSet<Likelihood>(Arrays.asList(this));
	}

	@Override
	public boolean isUsed() {
		return true;
	}// END: isUsed

	@Override
	public void setUsed() {
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

	
	
	
	
	
	
	
//	private LinkedList<Likelihood> getBranchLikelihoods() {
//
//		// linked list preserves order
//		LinkedList<Likelihood> loglikes = new LinkedList<Likelihood>();
//
//			for (NodeRef branch : treeModel.getNodes()) {
//
//				if (!treeModel.isRoot(branch)) {
//
//					int branchCategory = categoriesProvider.getBranchCategory(
//							treeModel, branch);
//					int index = (int) categoriesParameter
//							.getParameterValue(branchCategory);
//
//					Likelihood branchLikelihood = uniqueLikelihoods.get(index);
//
//					loglikes.add(branchLikelihood);
//					compoundModel.addModel(branchLikelihood.getModel());
//
//				}
//			}// END: branch loop
//
//		return loglikes;
//	}// END: getBranchLikelihoods
	
	
	
	
}// END: class
