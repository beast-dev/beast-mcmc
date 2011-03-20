package dr.evomodelxml.branchratemodel;

import dr.evomodel.branchratemodel.ArbitraryBranchRates;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.Parameter;
import dr.xml.*;

import java.util.logging.Logger;

/**
 */
public class ArbitraryBranchRatesParser extends AbstractXMLObjectParser {

    public static final String ARBITRARY_BRANCH_RATES = "arbitraryBranchRates";
    public static final String RATES = "rates";
    public static final String RECIPROCAL = "reciprocal";

    public String getParserName() {
        return ARBITRARY_BRANCH_RATES;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        TreeModel tree = (TreeModel) xo.getChild(TreeModel.class);

        XMLObject cxo = xo.getChild(RATES);

        Parameter rateCategoryParameter = (Parameter) cxo.getChild(Parameter.class);

        boolean reciprocal = xo.getAttribute(RECIPROCAL, false);

        final int numBranches = tree.getNodeCount() - 1;
        if (rateCategoryParameter.getDimension() != numBranches) {
            throw new XMLParseException("Invalid length for '" + rateCategoryParameter.getId() + "'\n" +
            "Should have length = " + numBranches);
        }

        Logger.getLogger("dr.evomodel").info("Using an scaled mixture of normals model.");
        Logger.getLogger("dr.evomodel").info("  rates = " + rateCategoryParameter.getDimension());
        Logger.getLogger("dr.evomodel").info("  NB: Make sure you have a prior on " + rateCategoryParameter.getId() + " and do not use this model in a treeLikelihood for sequence data");
        Logger.getLogger("dr.evomodel").info("  reciprocal = " + reciprocal);

        return new ArbitraryBranchRates(tree, rateCategoryParameter, reciprocal);
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "This element returns an arbitrary rate model." +
                "The branch rates are drawn from an arbitrary distribution determine by the prior.";
    }

    public Class getReturnType() {
        return ArbitraryBranchRates.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(TreeModel.class),
            new ElementRule(RATES, Parameter.class, "The rate parameter"),
            AttributeRule.newBooleanRule(RECIPROCAL, true),
    };


}
