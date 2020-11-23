package dr.evomodel.tree;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeUtils;
import dr.evolution.util.Taxa;
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
 * NOTE: current implementation requires paraphylyList to be disjoint set of tips that define paraphylies such that their union is the entire tree.
 * Could still implement option 2: keep track of branches already visited and traverse pre-order, but it may be more computationally expensive
 * todo: get statistic name right in logger
 * todo: add location option
 * todo: add optional weightByBranchTime
 */

public class ParaphylyRateStatistic extends TreeStatistic {
    public static final String PARAPHYLY_RATE_STATISTIC = "paraphylyRateStatistic";
    public static final String PARAPHYLY_LIST = "paraphylyList";
    private List<Taxa> paraphylySet;
    private DifferentiableBranchRates branchRateModel;
    private Tree tree;
    private Boolean weightByBranchTime;
    private int[] numNodesInClade;
    private List<NodeRef> MRCANodeList;
    private int dim;

    public ParaphylyRateStatistic(String name, DifferentiableBranchRates branchRateModel, List<Taxa> paraphylySet, Boolean weightByBranchTime, int dim) {
        super(name);
        this.branchRateModel = branchRateModel;
        this.paraphylySet = paraphylySet;
        this.tree = branchRateModel.getTree();
        this.weightByBranchTime = weightByBranchTime;
        this.dim = dim;
        this.MRCANodeList = new ArrayList<NodeRef>(dim);
        for (int i = 0; i < dim; i++) {
            MRCANodeList.add(null);
        }

        this.numNodesInClade = new int[dim];
        for (int i = 0; i < dim; i++) {
            numNodesInClade[i] = 2 * paraphylySet.get(i).getTaxonCount() - 1;
        }
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
        updateMRCAList();
        List<NodeRef> MRCANodeListComplement = new ArrayList<>(dim - 1); // all MRCA nodes except current paraphyly

        for (int j = 0; j < dim; j++) {
            if (j != i) {
                MRCANodeListComplement.add(MRCANodeList.get(j));
            }
        }

        double rateAverage = recurseToAccumulateRate(MRCANodeList.get(i), MRCANodeListComplement);
        return rateAverage / numNodesInClade[i];
    }

    private void updateMRCAList() {
        NodeRef node;
        for (int i = 0; i < dim; i++) {
            try {
                Set<String> leafSet = TreeUtils.getLeavesForTaxa(tree, paraphylySet.get(i));
                node = TreeUtils.getCommonAncestorNode(tree, leafSet);
                if (node == null) throw new RuntimeException("No clade found that contains " + leafSet);
                MRCANodeList.set(i, node);
            } catch (TreeUtils.MissingTaxonException e) {
                throw new RuntimeException("Missing taxon!");
            }
        }
    }

    private double recurseToAccumulateRate(NodeRef node, List<NodeRef> complement) {
        double total = 0.0;
        // curent default behavior includes "stem" of MRCA
        if (!tree.isExternal(node)) {
            if (!complement.contains(node)) {
                total += recurseToAccumulateRate(tree.getChild(node, 0), complement);
                total += recurseToAccumulateRate(tree.getChild(node, 1), complement);
            }
        }

        //ensures you don't add a root stem
        if (!complement.contains(node) && !tree.isRoot(node)) {
            total += branchRateModel.getUntransformedBranchRate(tree, node);
        }
        return total;
    }

    // **************************************************************
    // XMLObjectParser
    // **************************************************************

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return PARAPHYLY_RATE_STATISTIC;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            final String name = xo.getAttribute(Statistic.NAME, xo.getId());

            List<String> names = new ArrayList<>();


            DifferentiableBranchRates branchRateModel = (DifferentiableBranchRates) xo.getChild(DifferentiableBranchRates.class);

            List<Taxa> paraphylySet = new ArrayList<>();

            Taxa paraphyly;

            if (xo.hasChildNamed(PARAPHYLY_LIST)) {
                XMLObject cxo = xo.getChild(PARAPHYLY_LIST);
                for (int i = 0; i < cxo.getChildCount(); i++) {
                    paraphyly = (Taxa) cxo.getChild(i);
                    paraphylySet.add(paraphyly);
                }
            }

            int dim = paraphylySet.size();

            Boolean weightByBranchTime = false;

            ParaphylyRateStatistic paraphylyRateStatistic = new ParaphylyRateStatistic(name, branchRateModel, paraphylySet, weightByBranchTime, dim);
            return paraphylyRateStatistic;
        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
                new ElementRule(DifferentiableBranchRates.class),
                new ElementRule(PARAPHYLY_LIST, new XMLSyntaxRule[]{
                        new ElementRule(Taxa.class, 1, Integer.MAX_VALUE),
                }),
        };

        public String getParserDescription() {
            return "Reports weighted or unweighted average branch-rate parameter across complete disjoint set of user-specified paraphylies";
        }

        public Class getReturnType() {
            return ParaphylyRateStatistic.class;
        }
    };
}
