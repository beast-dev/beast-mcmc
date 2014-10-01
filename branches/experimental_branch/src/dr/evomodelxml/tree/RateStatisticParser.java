package dr.evomodelxml.tree;

import dr.evolution.tree.Tree;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.tree.RateStatistic;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.Statistic;
import dr.xml.*;

/**
 */
public class RateStatisticParser extends AbstractXMLObjectParser {

    public static final String RATE_STATISTIC = "rateStatistic";
    public static final String MODE = "mode";
    public static final String MEAN = "mean";
    public static final String VARIANCE = "variance";
    public static final String COEFFICIENT_OF_VARIATION = "coefficientOfVariation";

    public String getParserName() {
        return RATE_STATISTIC;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        final String name = xo.getAttribute(Statistic.NAME, xo.getId());
        final Tree tree = (Tree) xo.getChild(Tree.class);
        final BranchRateModel branchRateModel = (BranchRateModel) xo.getChild(BranchRateModel.class);

        final boolean internal = xo.getBooleanAttribute("internal");
        final boolean external = xo.getBooleanAttribute("external");

        if (!(internal || external)) {
            throw new XMLParseException("At least one of internal and external must be true!");
        }

        final String mode = xo.getStringAttribute(MODE);

        return new RateStatistic(name, tree, branchRateModel, external, internal, mode);
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "A statistic that returns the average of the branch rates";
    }

    public Class getReturnType() {
        return RateStatistic.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(TreeModel.class),
            new ElementRule(BranchRateModel.class),
            AttributeRule.newBooleanRule("internal"),
            AttributeRule.newBooleanRule("external"),
            new StringAttributeRule("mode", "This attribute determines how the rates are summarized, can be one of (mean, variance, coefficientOfVariance)", new String[]{MEAN, VARIANCE, COEFFICIENT_OF_VARIATION}, false),
            new StringAttributeRule("name", "A name for this statistic primarily for the purposes of logging", true),
    };

}
