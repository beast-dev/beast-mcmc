package dr.app.beagle.evomodel.branchmodel.lineagespecific;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeTrait;
import dr.evolution.tree.TreeTraitProvider;
import dr.inference.model.CompoundParameter;

/**
 * @author Filip Bielejec
 * @version $Id$
 * 
 */
public class BranchSpecificTrait implements TreeTraitProvider {

//	private CompoundParameter parameter;
//	private BranchSpecific branchSpecific;
	private Helper helper;
	
	public BranchSpecificTrait(final LineageSpecificBranchModel branchSpecific, final CompoundParameter parameter ) {
		
//		this.parameter = parameter;
//		this.branchSpecific = branchSpecific;
		helper = new Helper();
		
		TreeTrait<Double> uTrait = new TreeTrait.D() {

			@Override
			public String getTraitName() {
				return parameter.getId();
			}

			@Override
			public dr.evolution.tree.TreeTrait.Intent getIntent() {
				return Intent.BRANCH;
			}

			@Override
			public Double getTrait(Tree tree, NodeRef branch) {

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
