package dr.evomodel.treedatalikelihood.discrete.beastBasedDiscreteTreeLikelihood;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.util.Taxon;
import dr.evomodel.treedatalikelihood.DataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.xml.Reportable;

import java.util.Arrays;

public final class DiscretePostOrderReport implements Reportable {

    private final TreeDataLikelihood treeDataLikelihood;
    private final PostOrderMessageProvider discreteDelegate;

    public DiscretePostOrderReport(TreeDataLikelihood treeDataLikelihood) {
        this.treeDataLikelihood = treeDataLikelihood;

        DataLikelihoodDelegate delegate = treeDataLikelihood.getDataLikelihoodDelegate();
        if (delegate instanceof PostOrderMessageProvider) {
            this.discreteDelegate = (PostOrderMessageProvider) delegate;
        } else {
            throw new IllegalArgumentException("TreeDataLikelihood delegate does not provide post-order messages");
        }
    }

    @Override
    public String getReport() {

        treeDataLikelihood.getLogLikelihood();
        Tree tree = treeDataLikelihood.getTree();

        int categoryCount = discreteDelegate.getCategoryCount();
        int patternCount = discreteDelegate.getPatternCount();
        int stateCount = discreteDelegate.getStateCount();

        double[] start = new double[stateCount];
        double[] end = new double[stateCount];

        StringBuilder sb = new StringBuilder();
        sb.append("Discrete post-order partials\n");
        sb.append("  treeDataLikelihood = ").append(treeDataLikelihood.getId()).append('\n');
        sb.append("  nodes = ").append(tree.getNodeCount()).append('\n');
        sb.append("  patterns = ").append(patternCount).append('\n');
        sb.append("  categories = ").append(categoryCount).append('\n');
        sb.append("  states = ").append(stateCount).append('\n');

        for (int i = 0; i < tree.getNodeCount(); i++) {
            NodeRef node = tree.getNode(i);
            int nodeNumber = node.getNumber();
            sb.append("node ").append(nodeNumber);
            if (tree.isRoot(node)) {
                sb.append(" root");
            }
            if (tree.isExternal(node)) {
                Taxon taxon = tree.getNodeTaxon(node);
                sb.append(" taxon=").append(taxon == null ? "null" : taxon.getId());
            }
            sb.append('\n');

            for (int c = 0; c < categoryCount; c++) {
                for (int p = 0; p < patternCount; p++) {
                    discreteDelegate.getPostOrderBranchTopStandardInto(nodeNumber, c, p, start);
                    discreteDelegate.getPostOrderBranchBottomInto(nodeNumber, c, p, end);

                    sb.append("  category ").append(c)
                            .append(" pattern ").append(p)
                            .append(" start=").append(Arrays.toString(start))
                            .append(" end=").append(Arrays.toString(end))
                            .append('\n');
                }
            }
        }

        return sb.toString();
    }
}
