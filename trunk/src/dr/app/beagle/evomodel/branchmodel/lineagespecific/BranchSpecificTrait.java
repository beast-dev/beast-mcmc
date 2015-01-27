package dr.app.beagle.evomodel.branchmodel.lineagespecific;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.text.NumberFormat;

import dr.app.beagle.evomodel.branchmodel.BranchModel;
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

				int[] uCats = branchModel.getBranchModelMapping(branch).getOrder();
				int branchParameterIndex = uCats[0];
				
				//TODO: get the right parameter from the model, this assumes they are in the same order as models in Mapping
				double value = (Double) parameter.getParameter(branchParameterIndex).getParameterValue(0);
	
				System.out.println("FUBAR");
				
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
	
    //TODO: return annotated tree
	public String toString() {
		
		TreeLogger treeLogger = null;
		
		try {
			
//			File file = new File("annotated.tree");
//			file.deleteOnExit();
//			FileOutputStream out;
//			out = new FileOutputStream(file);

			
			 treeLogger = new TreeLogger(treeModel, //
					 null, //
                    null, //
                    new  TreeTraitProvider[] {this}, //
                    new TabDelimitedFormatter(System.out), // 
                    1, //
                    true, //
                     true, //
                     false, //
                     null, //
                     null);
					 
//					 new TreeLogger(treeModel, new TabDelimitedFormatter(System.out), 1, true, true, false);

			
			treeLogger.startLogging();
//			treeLogger.stopLogging();
			
		} catch (Exception e) {
			e.printStackTrace();
		}//END: try-catch

		return treeLogger.getTree().toString();

	}//END: toString
    
}//END: class
