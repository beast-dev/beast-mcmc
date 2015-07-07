package dr.app.beagle.evomodel.branchmodel.lineagespecific;

import java.util.LinkedList;
import java.util.List;

import dr.app.beagle.evomodel.branchmodel.BranchModel;
import dr.app.beagle.evomodel.sitemodel.SiteRateModel;
import dr.app.beagle.evomodel.substmodel.FrequencyModel;
import dr.app.beagle.evomodel.substmodel.SubstitutionModel;
import dr.app.beagle.evomodel.treelikelihood.BeagleTreeLikelihood;
import dr.evolution.alignment.PatternList;
import dr.evolution.tree.NodeRef;
import dr.evomodel.branchratemodel.BranchRateModel;
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
