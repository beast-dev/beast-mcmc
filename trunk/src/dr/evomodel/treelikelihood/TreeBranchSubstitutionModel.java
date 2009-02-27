package dr.evomodel.treelikelihood;

import dr.inference.model.Model;
import dr.inference.model.AbstractModel;
import dr.inference.model.Parameter;
import dr.evomodel.sitemodel.SiteModel;
import dr.evomodel.substmodel.SubstitutionModel;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evolution.tree.Tree;
import dr.evolution.tree.NodeRef;

/**
 * TreeBranchSubstitutionModel - provides a model for making substitution model epoches.
 * @author Marc A. Suchard
 */

public class TreeBranchSubstitutionModel extends AbstractModel {

    public TreeBranchSubstitutionModel(String name) {
        super(name);
    }

    public TreeBranchSubstitutionModel(String name, SiteModel siteModel, SubstitutionModel substModel, BranchRateModel branchModel) {
        super(name);
        this.siteModel = siteModel;
        this.substModel = substModel;
        this.branchModel = branchModel;

        if (siteModel != null)
            addModel(siteModel);
        if (substModel != null)
            addModel(substModel);
        if (branchModel != null)
            addModel(branchModel);
    }

    public void getTransitionProbabilities(Tree tree, NodeRef node, int rateCategory, double[] probs) {

         NodeRef parent = tree.getParent(node);

         final double branchRate = branchModel.getBranchRate(tree, node);

         // Get the operational time of the branch
         final double branchTime = branchRate * (tree.getNodeHeight(parent) - tree.getNodeHeight(node));

         if (branchTime < 0.0) {
             throw new RuntimeException("Negative branch length: " + branchTime);
         }

         double branchLength = siteModel.getRateForCategory(rateCategory) * branchTime;
         substModel.getTransitionProbabilities(branchLength, probs);
     } // getTransitionProbabilities

    protected void handleModelChangedEvent(Model model, Object object, int index) {
    }

    protected void handleParameterChangedEvent(Parameter parameter, int index, Parameter.ChangeType type) {
    }

    protected void storeState() {
    }

    protected void restoreState() {
    }

    protected void acceptState() {
    }

    protected SiteModel siteModel;
    protected SubstitutionModel substModel;
    protected BranchRateModel branchModel;

}
