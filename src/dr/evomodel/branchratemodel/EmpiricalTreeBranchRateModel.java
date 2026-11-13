package dr.evomodel.branchratemodel;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evomodel.tree.EmpiricalTreeDistributionModel;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;

/**
 * @author Marc A Suchard
 * @author Nidia Trovao
 */

public class EmpiricalTreeBranchRateModel extends AbstractBranchRateModel implements BranchRateModel {

    private final EmpiricalTreeDistributionModel trees;
    private final String rateAttributeName;

    private static final boolean CACHE = false;
    private final double[][] rateCache;

    public EmpiricalTreeBranchRateModel(String name,
                                        EmpiricalTreeDistributionModel trees,
                                        String rateAttributeName) {
        super(name);
        this.trees = trees;
        this.rateAttributeName = rateAttributeName;

        addModel(trees);

        if (CACHE) {
            rateCache = new double[trees.getTrees().size()][];
        } else {
            rateCache = null;
        }
    }

    @Override @SuppressWarnings("deprecated")
    public double getBranchRate(Tree tree, NodeRef node) {

        if (CACHE) {
            int currentIndex = trees.getCurrentTreeIndex();
            if (rateCache[currentIndex] == null) {
                rateCache[currentIndex] = getRatesAcrossTree(currentIndex);
            }
            return rateCache[currentIndex][node.getNumber()];
        } else {
            return extractRate(trees, node);
        }
    }

    @SuppressWarnings("deprecation")
    private double extractRate(Tree tree, NodeRef node) {
        return (double) tree.getNodeAttribute(node, rateAttributeName);
    }

    private double[] getRatesAcrossTree(int index) {

        Tree tree = trees.getTrees().get(index);
        NodeRef root = tree.getRoot();
        
        int length = tree.getNodeCount();
        double[] rates = new double[length];

        for (int i = 0; i < length; ++i) {
            NodeRef node = tree.getNode(i);
            if (node != root) {
                rates[node.getNumber()] = extractRate(tree, node);
            }
        }

        return rates;
    }

    @Override
    protected void handleModelChangedEvent(Model model, Object object, int index) {

    }

    @Override
    protected void storeState() {

    }

    @Override
    protected void restoreState() {

    }

    @Override
    protected void acceptState() {

    }

    @Override
    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {

    }
}
