package dr.app.beagle.evomodel.branchmodel.lineagespecific;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.text.NumberFormat;
import java.util.HashMap;
import java.util.Map;

import dr.app.beagle.evomodel.branchmodel.BranchModel;
import dr.app.bss.Utils;
import dr.evolution.tree.BranchRates;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeAttributeProvider;
import dr.evolution.tree.TreeTrait;
import dr.evolution.tree.TreeTraitProvider;
import dr.evomodel.branchratemodel.DefaultBranchRateModel;
import dr.evomodel.tree.TreeLogger;
import dr.evomodel.tree.TreeModel;
import dr.inference.loggers.LogFormatter;
import dr.inference.loggers.TabDelimitedFormatter;
import dr.inference.model.CompoundParameter;

/**
 * @author Filip Bielejec
 * @version $Id$
 * 
 */
public class BranchSpecificTrait implements TreeTraitProvider {

	private Helper helper;
	private TreeModel treeModel;
	
	public BranchSpecificTrait(
			TreeModel treeModel,
			final BranchModel branchModel,
			final CompoundParameter parameter ) {
		
		this.treeModel = treeModel;
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

				double value = 0.0;
				
				int[] uCats = branchModel.getBranchModelMapping(branch).getOrder();
				int category = uCats[0];

				//TODO: write some mechanism to get the right parameter from the substitution model
				value = (double) branchModel.getSubstitutionModels().get(category).getVariable(0).getValue(0);
				
//				System.out.println( 				branchModel.getSubstitutionModels().get(category).getVariable(0));
				System.out.println(value);
				
				
				System.out.println("--------------");
				
				
				
	
				return value;

			}//END: getTrait
		};
		
		helper.addTrait(uTrait);
		
	}//END: Constructor

    public TreeTrait[] getTreeTraits() {
        return helper.getTreeTraits();
    }

    public TreeTrait getTreeTrait(String key) {
        return helper.getTreeTrait(key);
    }

	// TODO: return annotated tree
	public String toString() {

		String annotatedTree = Tree.Utils.newick(treeModel, 
				new TreeTraitProvider[] { this });

		return annotatedTree;

	}// END: toString
    
}//END: class
