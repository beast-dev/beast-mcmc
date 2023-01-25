package dr.evomodelxml.branchratemodel;

import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.branchratemodel.BranchSpecificFixedEffects;
import dr.evomodel.branchratemodel.LocationScaledBranchRateModel;
import dr.evomodel.tree.TreeModel;
import dr.xml.*;

public class LocationScaledBranchRateModelParser extends AbstractXMLObjectParser {
    private static final String LOCATION_SCALED_BRANCH_RATE_MODEL = "locationScaledBranchRateModel";
    public String getParserName() {
        return LOCATION_SCALED_BRANCH_RATE_MODEL;
    }


    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        TreeModel tree = (TreeModel) xo.getChild(TreeModel.class);

        BranchRateModel branchRateModel = (BranchRateModel) xo.getChild(BranchRateModel.class);

        BranchSpecificFixedEffects location = (BranchSpecificFixedEffects) xo.getChild(BranchSpecificFixedEffects.class);

        return new LocationScaledBranchRateModel(tree, branchRateModel, location);
    }


    @Override
    public Class getReturnType() {
        return LocationScaledBranchRateModel.class;
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(TreeModel.class),
            new ElementRule(BranchRateModel.class),
            new ElementRule(BranchSpecificFixedEffects.class, "The location parameter"),
    };

    public String getParserDescription() {
        return "Returns a location scaled branch-rate model." +
                "All branch-rates are multiplied by the location parameter.";
    }

}

