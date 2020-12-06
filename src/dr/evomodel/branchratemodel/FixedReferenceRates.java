package dr.evomodel.branchratemodel;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.util.Taxon;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.MaskFromTree;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;
import dr.xml.*;

//import java.util.List;
import java.util.function.DoubleBinaryOperator;

/**
 * @author Alexander Fisher
 */

public class FixedReferenceRates extends AbstractBranchRateModel implements DifferentiableBranchRates{
    public static final String FIXED_REFERENCE_RATES = "fixedReferenceRates";
    public static final String FIXED_LENGTH = "fixedLength";

    private final TreeModel treeModel;
    private final BranchRateModel branchRateModel;
    private final DifferentiableBranchRates differentiableBranchRateModel;
    private final Taxon referenceTaxon;
    private final int fixedLength;
    private final MaskFromTree maskParameter;
//    private final Parameter rateParameter;
//    private List<NodeRef> nodeList;
    private NodeRef oneNode;
    private ArbitraryBranchRates arbitraryBranchRateModel;


    public FixedReferenceRates(String name, TreeModel treeModel, BranchRateModel branchRateModel, Taxon referenceTaxon, int fixedLength, MaskFromTree maskParameter) {
        super(name);
        this.treeModel = treeModel;
        this.branchRateModel = branchRateModel;
        this.referenceTaxon = referenceTaxon; //currently redundant
        this.fixedLength = fixedLength; //todo: add implementation for this optional parameter
        this.differentiableBranchRateModel = (branchRateModel instanceof DifferentiableBranchRates) ?
                (DifferentiableBranchRates) branchRateModel : null;
        this.maskParameter = maskParameter;

        checkDifferentiability(); //makes all the rest redundant
//        this.rateParameter = differentiableBranchRateModel.getRateParameter();
//        updateNodeList(treeModel, this.referenceTaxon);

        NodeRef node = updateNodeList(treeModel, this.referenceTaxon);
        assert(branchRateModel instanceof ArbitraryBranchRates);

        this.arbitraryBranchRateModel = (branchRateModel instanceof ArbitraryBranchRates) ?
                (ArbitraryBranchRates) branchRateModel : null;

//        arbitraryBranchRateModel.setBranchRate(treeModel, node, 1.0);


        addModel(treeModel);
        addModel(branchRateModel);

    }

    public double getUntransformedBranchRate(final Tree tree, final NodeRef node) {
//        if (nodeList.contains(node)) {
        updateNodeList(tree, referenceTaxon);
        if (node == oneNode) {
            return 1.0;
        }
        else {
            return arbitraryBranchRateModel.getUntransformedBranchRate(tree, node);
//            return 1.0;
        }
    }

    NodeRef updateNodeList(Tree tree, Taxon taxon) {

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
            node= tree.getParent(node);
        }

        oneNode = node;
        return oneNode;
//        this.nodeList.add(node);

//        int maskIndex;
//        for (int i = 0; i < nodeCache.length; i++) {
//            maskIndex = branchRates.getParameterIndexFromNodeNumber(nodeCache[i].getNumber());
//            this.setParameterValue(maskIndex, 0.0);
//        }
    }

//    private void setRateParameter(){
//        for(int i = 0; i < rateParameter.getDimension(); i++) {
//            if(maskParameter.getValue(i) == 0.0) {
//                rateParameter.setParameterValue(i, 1.0);
//            }
//        }
//    }

    public Tree getTree() {
        return treeModel;
    }

    @Override
    public double getBranchRateDifferential(final Tree tree, final NodeRef node) {
        checkDifferentiability();
        return differentiableBranchRateModel.getBranchRateDifferential(tree, node);
    }

    @Override
    public double getBranchRateSecondDifferential(Tree tree, NodeRef node) {
        checkDifferentiability();
        return differentiableBranchRateModel.getBranchRateSecondDifferential(tree, node);
    }

    @Override
    public Parameter getRateParameter() {
        checkDifferentiability();
        return differentiableBranchRateModel.getRateParameter();
    }

    @Override
    public int getParameterIndexFromNode(NodeRef node) {
        checkDifferentiability();
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
        checkDifferentiability();
        return differentiableBranchRateModel.mapReduceOverRates(map, reduce, initial);    }

    @Override
    public void forEachOverRates(NodeRateMap map) {
        checkDifferentiability();
        differentiableBranchRateModel.forEachOverRates(map);
    }

    @Override
    public double getPriorRateAsIncrement(Tree tree) {
        return 0;
    }

    @Override
    public double getBranchRate(Tree tree, NodeRef node) {
//        this.nodeList.clear();
//        updateNodeList(tree, this.referenceTaxon);

//        if(nodeList.contains(node)){
//        if (oneNode == node) {
//            return 1.0;
//        }
//        else {
//            return branchRateModel.getBranchRate(tree, node);
        return getUntransformedBranchRate(tree, node);
//        }
    }

    @Override
    protected void handleModelChangedEvent(Model model, Object object, int index) {
//        setRateParameter();
//        NodeRef node = updateNodeList(treeModel, this.referenceTaxon);
//        arbitraryBranchRateModel.setBranchRate(treeModel, node, 1.0);
        fireModelChanged();
    }

    @Override
    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
//        setRateParameter();
        fireModelChanged();
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

            MaskFromTree maskParameter = (MaskFromTree) xo.getChild(MaskFromTree.class);

            int fixedLength = xo.getAttribute(FIXED_LENGTH, 0);

//            int numBranches = tree.getNodeCount() - 1;

            FixedReferenceRates fixedReferenceRates = new FixedReferenceRates(FIXED_REFERENCE_RATES, tree, model, referenceTaxon, fixedLength, maskParameter);
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
                new ElementRule(MaskFromTree.class),
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
