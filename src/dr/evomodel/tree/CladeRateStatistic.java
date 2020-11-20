package dr.evomodel.tree;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeUtils;
import dr.evolution.util.Taxa;
import dr.evomodel.branchratemodel.ArbitraryBranchRates;
import dr.inference.model.Statistic;
import dr.xml.*;

import java.util.Set;

/**
 * @author Alexander Fisher
 */

/**
 * This statistic reports weighted and unweighted averages of arbitraryBranchRates across a user-specified clade.
 * NOTE: could potentially be refactored into CladeMeanAttributeStatistic?
 */
public class CladeRateStatistic extends TreeStatistic {
    public static final String CLADE_RATE_STATISTIC = "cladeRateStatistic";
    private Taxa subTreeLeafSet;
    private ArbitraryBranchRates branchRateModel;
    private Tree tree;
    private Boolean weightByBranchTime;
    private int numNodesInClade;
    private NodeRef MRCANode;

    // todo: currently BROKEN -- two options, 1. create MRCA set and traverse post order until you hit a different MRCA
    // todo: option 2: keep track of branches already visited and traverse pre-order
    public CladeRateStatistic(String name, ArbitraryBranchRates branchRateModel, Taxa subTreeLeafSet, Boolean weightByBranchTime) {
        super(name);
        this.branchRateModel = branchRateModel;
        this.subTreeLeafSet = subTreeLeafSet;
        this.tree = branchRateModel.getTree();
        this.weightByBranchTime = weightByBranchTime;
        this.numNodesInClade = 2 * subTreeLeafSet.getTaxonCount() - 2;
    }

    public void setTree(Tree tree) {
        this.tree = tree;
    }

    public Tree getTree() {
        return tree;
    }

    public int getDimension() {
        return 1;
    }

    @Override
    public double getStatisticValue(int dim) {
        NodeRef node;
        try {
            Set<String> leafSet = TreeUtils.getLeavesForTaxa(tree, subTreeLeafSet);
            node = TreeUtils.getCommonAncestorNode(tree, leafSet);
            if (node == null) throw new RuntimeException("No clade found that contains " + leafSet);
            MRCANode = node;
        } catch (TreeUtils.MissingTaxonException e) {
            throw new RuntimeException("Missing taxon!");
        }

        double rateAverage = recurseToAccumulateRate(node);

        return rateAverage / numNodesInClade;
    }

    private double recurseToAccumulateRate(NodeRef node) {
        double total = 0.0;

        // curent default behavior does not include stem
        // todo: set stem inclusion as optional attribute
        if (!tree.isExternal(node)) {
            total += recurseToAccumulateRate(tree.getChild(node, 0));
            total += recurseToAccumulateRate(tree.getChild(node, 1));
        }

        if (!tree.isRoot(node) && (node != MRCANode)) {
            total += branchRateModel.getUntransformedBranchRate(tree, node);
        }
        return total;
    }

    // **************************************************************
    // XMLObjectParser
    // **************************************************************

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return CLADE_RATE_STATISTIC;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            final String name = xo.getAttribute(Statistic.NAME, xo.getId());

            ArbitraryBranchRates branchRateModel = (ArbitraryBranchRates) xo.getChild(ArbitraryBranchRates.class);

            Taxa subTreeLeafSet = (Taxa) xo.getChild(Taxa.class);

            Boolean weightByBranchTime = false;

            // TODO: add optional weightByBranchTime
            CladeRateStatistic cladeRateStatistic = new CladeRateStatistic(name, branchRateModel, subTreeLeafSet, weightByBranchTime);
            return cladeRateStatistic;
        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
                new ElementRule(ArbitraryBranchRates.class),
                new ElementRule(Taxa.class),
        };

        public String getParserDescription() {
            return "Masks ancestral (off-root) branch to a specific reference taxon on a random tree";
        }

        public Class getReturnType() {
            return CladeRateStatistic.class;
        }
    };
}
