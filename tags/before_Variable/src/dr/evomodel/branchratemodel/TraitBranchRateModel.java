package dr.evomodel.branchratemodel;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.inference.model.Parameter;
import dr.inference.model.AbstractModel;
import dr.inference.model.Model;
import dr.xml.*;
import dr.evomodel.tree.TreeModel;

import java.util.logging.Logger;

/**
 * Takes the log rates at each node provided by a specified rate and give the branch rate as the average.
 *
 * @author Andrew Rambaut
 */
public class TraitBranchRateModel extends AbstractModel implements BranchRateModel {


    public static final String TRAIT_BRANCH_RATES = "traitBranchRates";
    public static final String TRAIT = "trait";
    public static final String RATE = "rate";
    public static final String RATIO = "ratio";

    private final String trait;
    private final Parameter rateParameter;
    private final Parameter ratioParameter;

    public TraitBranchRateModel(String trait, Parameter rateParameter, Parameter ratioParameter) {
        super(TRAIT_BRANCH_RATES);

        this.trait = trait;
        this.rateParameter = rateParameter;
        this.ratioParameter = ratioParameter;

        if (rateParameter != null) {
            addParameter(rateParameter);
        }

        if (ratioParameter != null) {
            addParameter(ratioParameter);
        }
    }

    public void handleModelChangedEvent(Model model, Object object, int index) {
    }

    protected final void handleParameterChangedEvent(Parameter parameter, int index, Parameter.ChangeType type) {
        fireModelChanged();
    }

    protected void storeState() {
    }

    protected void restoreState() {
    }

    protected void acceptState() {
    }


    public double getBranchRate(Tree tree, NodeRef node) {
        NodeRef parent = tree.getParent(node);
        if (parent == null) {
            throw new IllegalArgumentException("Root does not have a valid rate");
        }

        double scale = 1.0;
        double ratio = 1.0;

        if (rateParameter != null) {
            scale = rateParameter.getParameterValue(0);
        }

        if (ratioParameter != null) {
            ratio = ratioParameter.getParameterValue(0);
        }

        TreeModel treeModel = (TreeModel)tree;

        // get the log rate for the node and its parent
        double rate1 = ratio * treeModel.getMultivariateNodeTrait(node, trait)[0];
        double rate2 = ratio * treeModel.getMultivariateNodeTrait(parent, trait)[0];

        if (rate1 == rate2) {
            return scale * Math.exp(rate1);
        }

        double rate = scale * (Math.exp(rate2) - Math.exp(rate1)) / (rate2 - rate1);

        return rate;
    }

    public String getBranchAttributeLabel() {
        return "rate";
    }

    public String getAttributeForBranch(Tree tree, NodeRef node) {
        return Double.toString(getBranchRate(tree, node));
    }

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return TRAIT_BRANCH_RATES;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            String trait = xo.getStringAttribute(TRAIT);

            Logger.getLogger("dr.evomodel").info("Using trait, " + trait + ", as log rate estimates.");

            Parameter rateParameter = (Parameter) xo.getElementFirstChild(RATE);
            Parameter ratioParameter = (Parameter) xo.getElementFirstChild(RATIO);

            return new TraitBranchRateModel(trait, rateParameter, ratioParameter);
        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public String getParserDescription() {
            return
                    "This element returns an trait rate model." +
                            "The branch rates are an average of the rates provided by a node trait.";
        }

        public Class getReturnType() {
            return TraitBranchRateModel.class;
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
	            AttributeRule.newStringRule(TRAIT, false, "The name of the trait that provides the log rates at nodes"),
                new ElementRule(RATE, Parameter.class, "The rate parameter", true),
                new ElementRule(RATIO, Parameter.class, "The ratio parameter", true),
        };
    };


}