package dr.evomodel.continuous;

import dr.evolution.tree.BranchAttributeProvider;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evomodel.tree.TreeModel;

/**
 * @author Marc Suchard
 */
public abstract class BivariateTraitBranchAttributeProvider implements BranchAttributeProvider {

    public static final String FORMAT = "%5.4f";

    public BivariateTraitBranchAttributeProvider(SampledMultivariateTraitLikelihood traitLikelihood) {

        traitName = traitLikelihood.getTraitName();
        label = traitName + extensionName();
        treeModel = traitLikelihood.getTreeModel();

        double[] rootTrait = treeModel.getMultivariateNodeTrait(treeModel.getRoot(), traitName);
        if (rootTrait.length != 2)
            throw new RuntimeException("BivariateTraitBranchAttributeProvider only works for 2D traits");
    }

    protected abstract String extensionName();

    protected double branchFunction(double[] startValue, double[] endValue, double startTime, double endTime) {
         return convert(endValue[0]-startValue[0], endValue[1] - startValue[1], startTime - endTime);
    }

    protected abstract double convert(double latDifference, double longDifference, double timeDifference);

    public String getBranchAttributeLabel() {
        return label;
    }

    public String getAttributeForBranch(Tree tree, NodeRef node) {

        if (tree != treeModel)
            throw new RuntimeException("Bad bug.");

        NodeRef parent = tree.getParent(node);
        double[] startTrait = treeModel.getMultivariateNodeTrait(parent, traitName);
        double[] endTrait = treeModel.getMultivariateNodeTrait(node, traitName);
        double startTime = tree.getNodeHeight(parent);
        double endTime = tree.getNodeHeight(node);

        return String.format(BivariateTraitBranchAttributeProvider.FORMAT,
                branchFunction(startTrait, endTrait, startTime, endTime));

    }

    protected TreeModel treeModel;
    protected String traitName;
    protected String label;

}
