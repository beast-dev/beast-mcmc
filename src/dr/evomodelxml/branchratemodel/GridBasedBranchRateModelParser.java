package dr.evomodelxml.branchratemodel;

import dr.evomodel.branchratemodel.GridBasedBranchRateModel;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.Parameter;

import dr.xml.*;

public class GridBasedBranchRateModelParser extends AbstractXMLObjectParser {
    public final String PARSER_NAME = "gridBasedBranchRateModel";
    private final String GRID_POINTS = "gridPoints"; // TODO make this more flexible; create gridPoints class
    private final String RATES_FUNCTION = "levelSpecificRates";

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

//        int dim = xo.getIntegerAttribute(DIMENSION);

        TreeModel tree = (TreeModel) xo.getChild(TreeModel.class);
        Parameter levelSpecificRates = (Parameter) xo.getElementFirstChild(RATES_FUNCTION);

        Parameter gridPoints = (Parameter) xo.getElementFirstChild(GRID_POINTS);
        if (gridPoints.getDimension() != levelSpecificRates.getDimension() - 1) {
            throw new XMLParseException("The number of grid points should be one less than the number of rates");
        }

        return new GridBasedBranchRateModel(tree, gridPoints, levelSpecificRates);
    }

    public XMLSyntaxRule[] getSyntaxRules() { return rules; }

    private final XMLSyntaxRule[] rules = {
//            AttributeRule.newIntegerRule(DIMENSION),
            new ElementRule(TreeModel.class),
            new ElementRule(GRID_POINTS, Parameter.class),
            new ElementRule(RATES_FUNCTION, Parameter.class)
    };

    @Override
    public String getParserDescription() {
        return "Parser for grid based branch rate model";
    }

    @Override
    public Class getReturnType() {
        return GridBasedBranchRateModel.class;
    }

    @Override
    public String getParserName() {
        return PARSER_NAME;
    }
}
