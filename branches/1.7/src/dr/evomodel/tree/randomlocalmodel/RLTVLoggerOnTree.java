package dr.evomodel.tree.randomlocalmodel;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeTrait;

/**
 * @author Marc A. Suchard
 */
public class RLTVLoggerOnTree implements TreeTrait<Double> {

    public static final String TRAIT_NAME = "changed";

    public RLTVLoggerOnTree(RandomLocalTreeVariable treeVariable) {
        this.treeVariable = treeVariable;
    }

    public String getTraitName() {
        return TRAIT_NAME;
    }

    public Intent getIntent() {
        return Intent.BRANCH;
    }

    public Class getTraitClass() {
        return Double.class;
    }

    public boolean getLoggable() {
        return true;
    }

    public Double getTrait(final Tree tree, final NodeRef node) {
        return treeVariable.isVariableSelected(tree, node) ? 1.0 : 0.0;
    }

    public String getTraitString(final Tree tree, final NodeRef node) {
        return Double.toString(getTrait(tree, node));
    }

    private RandomLocalTreeVariable treeVariable;
}
