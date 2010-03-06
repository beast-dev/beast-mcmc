package dr.evomodel.continuous;

import dr.evolution.tree.BranchAttributeProvider;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;

/**
 * @author Marc Suchard
 */
public abstract class BivariateTraitBranchAttributeProvider implements BranchAttributeProvider {

    public static final String FORMAT = "%5.4f";

    public BivariateTraitBranchAttributeProvider(AbstractMultivariateTraitLikelihood traitLikelihood) {

        traitName = traitLikelihood.getTraitName();
        label = traitName + extensionName();

        double[] rootTrait = traitLikelihood.getRootNodeTrait();
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

        if (tree != traitLikelihood.getTreeModel())
            throw new RuntimeException("Bad bug.");

        NodeRef parent = tree.getParent(node);
        double[] startTrait = traitLikelihood.getTraitForNode(tree, parent, traitName);
        double[] endTrait = traitLikelihood.getTraitForNode(tree, node, traitName);
        double startTime = tree.getNodeHeight(parent);
        double endTime = tree.getNodeHeight(node);

        return String.format(BivariateTraitBranchAttributeProvider.FORMAT,
                branchFunction(startTrait, endTrait, startTime, endTime));

    }

    protected AbstractMultivariateTraitLikelihood traitLikelihood;
    protected String traitName;
    protected String label;

}
