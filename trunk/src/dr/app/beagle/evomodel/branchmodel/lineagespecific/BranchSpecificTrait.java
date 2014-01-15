package dr.app.beagle.evomodel.branchmodel.lineagespecific;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeTrait;
import dr.evolution.tree.TreeTraitProvider;
import dr.evolution.tree.TreeTrait.Intent;
import dr.evolution.tree.TreeTraitProvider.Helper;
import dr.inference.model.CompoundParameter;

public class BranchSpecificTrait implements  TreeTraitProvider {

	private CompoundParameter parameter;
	private BranchSpecific branchSpecific;
	private Helper helper;
	
	/*
	
	<fancyNewTreeTraitProvider id=“omegaPerBranch”>  <!— returns parameter(i) per branch —>
	<branchCategoryModel idref=“filipModel”/>  <!— returns 0, 1, 2, etc per branch —>
	<compoundParameter>
		<parameter idref=“omega1”/>
		<parameter idref=“omega2”/>   <!— implicit contract that parameters are listed in correct order,  could make explicit in XML at some point with an attribute —>
	</compoundParameter>
< fancyNewTreeTraitProvider>

<fancyNewTreeTraitProvider id=“kappaPerBranch”>  <!— returns parameter(i) per branch —>
	<branchCategoryModel idref=“filipModel”/>  <!— returns 0, 1, 2, etc per branch —>
	<compoundParameter>
		<parameter idref=“kappa1”/>
		<parameter idref=“kaapa2”/>   <!— implicit contract that parameters are listed in correct order,  could make explicit in XML at some point with an attribute —>
	</compoundParameter>
< fancyNewTreeTraitProvider>

<treeLog>  <!— automagically exports trees for FigTree or TreeAnnonator.
	<fancyNewTreeTraitProvider idref=“omegaPerBranch”/>
	<fancyNewTreeTraitProvider idref=“kappaPerBranch”/>
</treeLog>
	
	*/
	
	public BranchSpecificTrait(final BranchSpecific branchSpecific, final CompoundParameter parameter ) {
		
		this.parameter = parameter;
		this.branchSpecific = branchSpecific;
		helper = new Helper();
		
		TreeTrait<Double> uTrait = new TreeTrait.D() {

			@Override
			public String getTraitName() {
				
				
//				System.err.println(parameter.getParameterName());
				
				return parameter.getId();
			}

			@Override
			public dr.evolution.tree.TreeTrait.Intent getIntent() {
				return Intent.BRANCH;
			}

			@Override
			public Double getTrait(Tree tree, NodeRef branch) {

				branchSpecific.getBranchModelMapping(branch);
				
				int[] uCats = branchSpecific.getBranchModelMapping(branch).getOrder();
				int branchParameterIndex = uCats[0];
				
				//TODO: get the right parameter from the model, this assumes they are in the same order as models in Mapping
				double value = (Double) parameter.getParameter(branchParameterIndex).getParameterValue(0);
	
				return value;

			}
		};
		
		helper.addTrait(uTrait);
		
	}//END: Constructor

    public TreeTrait[] getTreeTraits() {
        return helper.getTreeTraits();
    }

    public TreeTrait getTreeTrait(String key) {
        return helper.getTreeTrait(key);
    }
	
}//END: class
