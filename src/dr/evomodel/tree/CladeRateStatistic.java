package dr.evomodel.tree;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeUtils;
import dr.evolution.util.Taxa;
import dr.evomodel.branchratemodel.AbstractBranchRateModel;
import dr.evomodel.branchratemodel.ArbitraryBranchRates;
import dr.evomodel.branchratemodel.DifferentiableBranchRates;
import dr.inference.model.Statistic;
import dr.xml.*;

import java.util.ArrayList;
import java.util.List;
import java.util.Set;

/**
 * @author Alexander Fisher
 */

/**
 * This statistic reports weighted and unweighted averages of arbitraryBranchRates across a all user-specified clades in a tree.
 * NOTE: current implementation requires cladeList to be disjoint set of clades such that their union is the entire tree
 * could still implement option 2: keep track of branches already visited and traverse pre-order, but it would be more computationally expensive
 * todo: get statistic name right in logger
 * todo: rename everything to "paraphyletic'
 */

public class CladeRateStatistic extends TreeStatistic {
    public static final String CLADE_RATE_STATISTIC = "cladeRateStatistic";
    public static final String CLADE_LIST = "cladeList";
    private List<Taxa> cladeSet;
    private DifferentiableBranchRates branchRateModel;
    private Tree tree;
    private Boolean weightByBranchTime;
    private int[] numNodesInClade;
    private List<NodeRef> MRCANodeList;
    private int dim;

    public CladeRateStatistic(String name, DifferentiableBranchRates branchRateModel, List<Taxa> cladeSet, Boolean weightByBranchTime, int dim) {
        super(name);
        this.branchRateModel = branchRateModel;
        this.cladeSet = cladeSet;
        this.tree = branchRateModel.getTree();
        this.weightByBranchTime = weightByBranchTime;
        this.dim = dim;
        //todo: update MRCA index list on every tree change
        this.MRCANodeList = new ArrayList<>(dim);
        this.numNodesInClade = new int[dim];
        for (int i = 0; i < dim; i++) {
            numNodesInClade[i] = 2 * cladeSet.get(i).getTaxonCount() - 2;
        }
        updateMRCAList();
    }

    public void setTree(Tree tree) {
        this.tree = tree;
    }

    public Tree getTree() {
        return tree;
    }

    public int getDimension() {
        return dim;
    }

    @Override
    public double getStatisticValue(int i) {
//        try {
//            Set<String> leafSet = TreeUtils.getLeavesForTaxa(tree, cladeSet.get(i));
//            node = TreeUtils.getCommonAncestorNode(tree, leafSet);
//            if (node == null) throw new RuntimeException("No clade found that contains " + leafSet);
//            MRCANode = node;
//        } catch (TreeUtils.MissingTaxonException e) {
//            throw new RuntimeException("Missing taxon!");
//        }

        double rateAverage = recurseToAccumulateRate(MRCANodeList.get(i));
        return rateAverage / numNodesInClade[i];
    }

    private void updateMRCAList() {
        NodeRef node;
        for (int i = 0; i < dim; i++) {
            try {
                Set<String> leafSet = TreeUtils.getLeavesForTaxa(tree, cladeSet.get(i));
                node = TreeUtils.getCommonAncestorNode(tree, leafSet);
                if (node == null) throw new RuntimeException("No clade found that contains " + leafSet);
                MRCANodeList.add(node);
            } catch (TreeUtils.MissingTaxonException e) {
                throw new RuntimeException("Missing taxon!");
            }
        }
    }

    private double recurseToAccumulateRate(NodeRef node) {
        double total = 0.0;

        // curent default behavior does not include stem
        // todo: set stem inclusion as optional attribute
        if (!tree.isExternal(node)) {
            total += recurseToAccumulateRate(tree.getChild(node, 0));
            total += recurseToAccumulateRate(tree.getChild(node, 1));
        }

        if (!tree.isRoot(node) && (!MRCANodeList.contains(node))) {
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

            List<String> names = new ArrayList<>();


            DifferentiableBranchRates branchRateModel = (DifferentiableBranchRates) xo.getChild(DifferentiableBranchRates.class);

            List<Taxa> cladeSet = new ArrayList<>();

            Taxa clade;

            if (xo.hasChildNamed(CLADE_LIST)) {
                XMLObject cxo = xo.getChild(CLADE_LIST);
                for (int i = 0; i < cxo.getChildCount(); i++) {
                    clade = (Taxa) cxo.getChild(i);
                    cladeSet.add(clade);
                }
            }

//            Taxa subTreeLeafSet = (Taxa) xo.getChild(Taxa.class);

            int dim = cladeSet.size();

            Boolean weightByBranchTime = false;

            // TODO: add optional weightByBranchTime
            CladeRateStatistic cladeRateStatistic = new CladeRateStatistic(name, branchRateModel, cladeSet, weightByBranchTime, dim);
            return cladeRateStatistic;
        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
                new ElementRule(DifferentiableBranchRates.class),
                new ElementRule(CLADE_LIST, new XMLSyntaxRule[]{
                        new ElementRule(Taxa.class, 1, Integer.MAX_VALUE),
                }),
        };

        public String getParserDescription() {
            return "Masks ancestral (off-root) branch to a specific reference taxon on a random tree";
        }

        public Class getReturnType() {
            return CladeRateStatistic.class;
        }
    };
}
