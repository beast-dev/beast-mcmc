package dr.evomodel.branchratemodel;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.util.Taxon;
import dr.evomodel.tree.TreeModel;
import dr.inference.loggers.Loggable;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;
import dr.xml.*;

//import java.util.List;
import java.util.function.DoubleBinaryOperator;

/**
 * @author Alexander Fisher
 */

public class FixedReferenceRates extends AbstractBranchRateModel implements DifferentiableBranchRates, Loggable {
    public static final String FIXED_REFERENCE_RATES = "fixedReferenceRates";
    public static final String FIXED_LENGTH = "fixedLength";

    private final TreeModel treeModel;
    private final DifferentiableBranchRates differentiableBranchRateModel;
    private final Taxon referenceTaxon;
    private final int fixedLength;
    //    private List<NodeRef> nodeList;
    private NodeRef oneNode; // node to be set to 1.0
    private boolean nodeKnown = false;
    private boolean storedNodeKnown;
    private NodeRef storedOneNode;

    public FixedReferenceRates(String name, TreeModel treeModel, BranchRateModel branchRateModel, Taxon referenceTaxon, int fixedLength) {
        super(name);
        this.treeModel = treeModel;
        this.referenceTaxon = referenceTaxon;
        this.fixedLength = fixedLength; //todo: add implementation for this optional parameter
        this.differentiableBranchRateModel = (branchRateModel instanceof DifferentiableBranchRates) ?
                (DifferentiableBranchRates) branchRateModel : null;

//        if (differentiableBranchRateModel instanceof ArbitraryBranchRates) {
//            ArbitraryBranchRates arbitraryBranchRateModel = (ArbitraryBranchRates) differentiableBranchRateModel;
//            if ((arbitraryBranchRateModel.getTransform() instanceof ArbitraryBranchRates.BranchRateTransform.None) == false) {
//                throw new RuntimeException("Transformed arbitrary branch rates not yet implemented.");
//            }
//        }

        checkDifferentiability();

        updateNodeList(treeModel, this.referenceTaxon);

        addModel(treeModel);
        addModel(branchRateModel);

    }

    public double getUntransformedBranchRate(final Tree tree, final NodeRef node) {

        if (!nodeKnown) {
            updateNodeList(tree, referenceTaxon);
        }
        if (node == oneNode) {
            return 1.0;
        } else {
            return differentiableBranchRateModel.getUntransformedBranchRate(tree, node);
        }
    }

    void updateNodeList(Tree tree, Taxon taxon) {
        nodeKnown = true;
        int nodeNumber = tree.getTaxonIndex(taxon.getId());
        NodeRef node = tree.getNode(nodeNumber);
        NodeRef root = tree.getRoot();
//        NodeRef[] nodeCache = new NodeRef[maskLength + 1];

        //set node cache to all be the same node
//        for (int i = 0; i < nodeCache.length; i++) {
//            nodeCache[i] = node;
//        }

//        int update = 0;
        while (tree.getParent(node) != root) {
//            update = update + 1;
//            for (int i = maskLength; i > 0; i--) {
//                if (update > i) {
//                    nodeCache[i] = nodeCache[i - 1];
//                }
//            }
            node = tree.getParent(node);
        }

        oneNode = node;
    }

    public Tree getTree() {
        return treeModel;
    }

    @Override
    public double getBranchRateDifferential(final Tree tree, final NodeRef node) {
        return differentiableBranchRateModel.getBranchRateDifferential(tree, node);
    }

    @Override
    public double getBranchRateSecondDifferential(Tree tree, NodeRef node) {
        return differentiableBranchRateModel.getBranchRateSecondDifferential(tree, node);
    }

    @Override
    public Parameter getRateParameter() {
        return differentiableBranchRateModel.getRateParameter();
    }

    @Override
    public int getParameterIndexFromNode(NodeRef node) {
        return differentiableBranchRateModel.getParameterIndexFromNode(node);
    }

    private void checkDifferentiability() {
        if (differentiableBranchRateModel == null) {
            throw new RuntimeException("Non-differentiable base BranchRateModel");
        }
    }

    @Override
    public ArbitraryBranchRates.BranchRateTransform getTransform() {
        throw new RuntimeException("Not yet implemented");
    }

    @Override
    public double[] updateGradientLogDensity(double[] gradient, double[] value, int from, int to) {
        return differentiableBranchRateModel.updateGradientLogDensity(gradient, value, from, to);
    }

    @Override
    public double[] updateDiagonalHessianLogDensity(double[] diagonalHessian, double[] gradient, double[] value, int from, int to) {
        throw new RuntimeException("Not yet implemented");
    }

    @Override
    public double mapReduceOverRates(NodeRateMap map, DoubleBinaryOperator reduce, double initial) {
        return differentiableBranchRateModel.mapReduceOverRates(map, reduce, initial);
    }

    @Override
    public void forEachOverRates(NodeRateMap map) {
        differentiableBranchRateModel.forEachOverRates(map);
    }

    @Override
    public double getPriorRateAsIncrement(Tree tree) {
        return 0;
    }

    @Override
    public double getBranchRate(Tree tree, NodeRef node) {
        if (!nodeKnown) {
            updateNodeList(tree, referenceTaxon);
        }
        if (node == oneNode) {
            return 1.0;
        } else {
            return differentiableBranchRateModel.getBranchRate(tree, node);
        }
    }

    @Override
    protected void handleModelChangedEvent(Model model, Object object, int index) {
        if (model == differentiableBranchRateModel) {
            fireModelChanged();
        } else if (model == treeModel) {
            nodeKnown = false;
        } else {
            throw new RuntimeException("Should only watch branchRates or treeModel");
        }
    }

    @Override
    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        throw new RuntimeException("Should not be variable changed event");
    }

    @Override
    protected void storeState() {
        storedNodeKnown = nodeKnown;
        storedOneNode = oneNode;
    }

    @Override
    protected void restoreState() {
        nodeKnown = storedNodeKnown;
        oneNode = storedOneNode;
    }

    @Override
    protected void acceptState() {

    }

//    @Override
//    LogColumn[] getColumns() {
//        return new LogColumn[differentiableBranchRateModel.getRateParameter().getDimension()] {
//            new LikelihoodColumn(getId())
//        }
//    }

    // **************************************************************
    // XMLObjectParser
    // **************************************************************

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return FIXED_REFERENCE_RATES;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            TreeModel tree = (TreeModel) xo.getChild(TreeModel.class);

            BranchRateModel model = (BranchRateModel) xo.getChild(BranchRateModel.class);

            Taxon referenceTaxon = (Taxon) xo.getChild(Taxon.class);

            int fixedLength = xo.getAttribute(FIXED_LENGTH, 0);

            FixedReferenceRates fixedReferenceRates = new FixedReferenceRates(FIXED_REFERENCE_RATES, tree, model, referenceTaxon, fixedLength);
            return fixedReferenceRates;
        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
                new ElementRule(TreeModel.class),
                new ElementRule(BranchRateModel.class),
                new ElementRule(Taxon.class),
                AttributeRule.newStringRule(FIXED_LENGTH, true),
        };

        public String getParserDescription() {
            return "Fixes ancestral off-root branch (and optional addnl branches) to 1 with reference to a user-specified taxon.";
        }

        public Class getReturnType() {
            return FixedReferenceRates.class;
        }
    };
}
