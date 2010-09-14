package dr.evomodelxml.branchratemodel;

import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.branchratemodel.CompoundBranchRateModel;
import dr.xml.*;

import java.util.ArrayList;
import java.util.List;
import java.util.logging.Logger;

/**
 */
public class CompoundBranchRateModelParser extends AbstractXMLObjectParser {

    public static final String COMPOUND_BRANCH_RATE_MODEL = "compoundBranchRateModel";

    public String getParserName() {
        return COMPOUND_BRANCH_RATE_MODEL;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        List<BranchRateModel> branchRateModels = new ArrayList<BranchRateModel>();
        for (int i = 0; i < xo.getChildCount(); i++) {
            if (xo.getChild(i) instanceof BranchRateModel) {
                branchRateModels.add((BranchRateModel) xo.getChild(i));
            } else {

                Object rogueElement = xo.getChild(i);

                throw new XMLParseException("An element (" + rogueElement
                        + ") which is not a branchRateModel has been added to a " + COMPOUND_BRANCH_RATE_MODEL + " element");
            }

        }

        Logger.getLogger("dr.evomodel").info("Creating a compound branch rate model of " + branchRateModels.size() + " sub-models");

        return new CompoundBranchRateModel(branchRateModels);
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "This element provides a strict clock model. " +
                "All branches have the same rate of molecular evolution.";
    }

    public Class getReturnType() {
        return CompoundBranchRateModel.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
            new ElementRule(BranchRateModel.class, "The component branch rate models", 1, Integer.MAX_VALUE),
    };

}
