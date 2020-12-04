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

    public FixedReferenceRates(TreeModel treeModel, BranchRateModel branchRateModel, Taxon referenceTaxon, int fixedLength) {
        super(this.FIXED_REFERENCE_RATES);
        this.treeModel = treeModel;
        this.branchRateModel = branchRateModel;
        this.referenceTaxon = referenceTaxon;
        this.fixedLength = fixedLength;
        this.differentiableBranchRateModel = (branchRateModel instanceof DifferentiableBranchRates) ?
                (DifferentiableBranchRates) branchRateModel : null;
    }

    @Override
    public double getBranchRateDifferential(Tree tree, NodeRef node) {
        return 0;
    }

    @Override
    public double getBranchRateSecondDifferential(Tree tree, NodeRef node) {
        return 0;
    }

    @Override
    public Parameter getRateParameter() {
        return null;
    }

    @Override
    public int getParameterIndexFromNode(NodeRef node) {
        return 0;
    }

    @Override
    public ArbitraryBranchRates.BranchRateTransform getTransform() {
        return null;
    }

    @Override
    public double[] updateGradientLogDensity(double[] gradient, double[] value, int from, int to) {
        return new double[0];
    }

    @Override
    public double[] updateDiagonalHessianLogDensity(double[] diagonalHessian, double[] gradient, double[] value, int from, int to) {
        return new double[0];
    }

    @Override
    public double mapReduceOverRates(NodeRateMap map, DoubleBinaryOperator reduce, double initial) {
        return 0;
    }

    @Override
    public void forEachOverRates(NodeRateMap map) {

    }

    @Override
    public double getPriorRateAsIncrement(Tree tree) {
        return 0;
    }

    @Override
    public double getBranchRate(Tree tree, NodeRef node) {
        return 0;
    }

    @Override
    protected void handleModelChangedEvent(Model model, Object object, int index) {

    }

    @Override
    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {

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

            int fixedLength = xo.getAttribute(FIXED_LENGTH, 0);

//            int numBranches = tree.getNodeCount() - 1;

            FixedReferenceRates fixedReferenceRates = new FixedReferenceRates(tree, model, referenceTaxon, fixedLength);
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
