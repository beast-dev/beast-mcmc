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
 */

public class ParaphylyRateStatistic extends TreeStatistic {
    public static final String PARAPHYLY_RATE_STATISTIC = "paraphylyRateStatistic";
    public static final String PARAPHYLY_LIST = "paraphylyList";
    public static final String WEIGHTING = "weighting";

    private List<Taxa> paraphylySet;
    private DifferentiableBranchRates branchRateModel;
    private Tree tree;
    private double totalTime;
    private List<NodeRef> MRCANodeList;
    private int dim;
    private BranchWeighting branchWeighting;

    public ParaphylyRateStatistic(String name, DifferentiableBranchRates branchRateModel, List<Taxa> paraphylySet, BranchWeighting branchWeighting, int dim) {
        super(name);
        this.branchRateModel = branchRateModel;
        this.paraphylySet = paraphylySet;
        this.tree = branchRateModel.getTree();
        this.branchWeighting = branchWeighting;
        this.dim = dim;
        this.MRCANodeList = new ArrayList<NodeRef>(dim);
        for (int i = 0; i < dim; i++) {
            MRCANodeList.add(null);
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

        this.totalTime = 0.0;
        double rateAverage = recurseToAccumulateRate(MRCANodeList.get(i), MRCANodeListComplement);
        return rateAverage / totalTime;
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
            total += branchWeighting.getBranchRate(branchRateModel, tree, node);
            this.totalTime += branchWeighting.getDenominator(tree, node);
        }
        return total;
    }

    public enum BranchWeighting {

        NONE("none") {
            @Override
            double getBranchRate(DifferentiableBranchRates branchRateModel, Tree tree, NodeRef node) {
                return branchRateModel.getUntransformedBranchRate(tree, node);
            }

            @Override
            double getDenominator(Tree tree, NodeRef node) {
                return 1.0;
            }
        },

        BY_TIME("byTime") {
            @Override
            double getBranchRate(DifferentiableBranchRates branchRateModel, Tree tree, NodeRef node) {
                return branchRateModel.getUntransformedBranchRate(tree, node) * tree.getBranchLength(node);
            }

            @Override
            double getDenominator(Tree tree, NodeRef node) {
                return tree.getBranchLength(node);
            }
        };

        BranchWeighting(String name) {
            this.name = name;
        }

        public String getName() {
            return name;
        }

        private final String name;

        abstract double getBranchRate(DifferentiableBranchRates branchRateModel, Tree tree, NodeRef node);

        abstract double getDenominator(Tree tree, NodeRef node);

        public static BranchWeighting parse(String name) {
            for (BranchWeighting weighting : BranchWeighting.values()) {
                if (weighting.getName().equalsIgnoreCase(name)) {
                    return weighting;
                }
            }
            return null;
        }
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

            BranchWeighting branchWeighting = parseWeighting(xo);

            ParaphylyRateStatistic paraphylyRateStatistic = new ParaphylyRateStatistic(name, branchRateModel, paraphylySet, branchWeighting, dim);
            return paraphylyRateStatistic;
        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        private ParaphylyRateStatistic.BranchWeighting parseWeighting(XMLObject xo)
                throws XMLParseException {

            String defaultName = ParaphylyRateStatistic.BranchWeighting.NONE.getName();
            String name = xo.getAttribute(WEIGHTING, defaultName);

            ParaphylyRateStatistic.BranchWeighting weighting =
                    ParaphylyRateStatistic.BranchWeighting.parse(name);

            if (weighting == null) {
                throw new XMLParseException("Unknown weighting type");
            }

            return weighting;
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
                new ElementRule(DifferentiableBranchRates.class),
                new ElementRule(PARAPHYLY_LIST, new XMLSyntaxRule[]{
                        new ElementRule(Taxa.class, 1, Integer.MAX_VALUE),
                }),
                AttributeRule.newStringRule(WEIGHTING, true),

        };

        public String getParserDescription() {
            return "Reports weighted or unweighted average branch-rate parameter across complete disjoint set of user-specified paraphylies";
        }

        public Class getReturnType() {
            return ParaphylyRateStatistic.class;
        }
    };
}
