package dr.evomodelxml.coalescent;

import dr.evolution.util.Units;
import dr.evomodel.coalescent.PiecewisePopulationModel;
import dr.evomodel.coalescent.ScaledPiecewiseModel;
import dr.evomodel.tree.TreeModel;
import dr.evoxml.util.XMLUnits;
import dr.inference.model.Parameter;
import dr.xml.*;

/**
 * Parses an element from an DOM document into a PiecewisePopulation.
 */
public class ScaledPiecewiseModelParser extends AbstractXMLObjectParser {

    public static final String PIECEWISE_POPULATION = "scaledPiecewisePopulation";
    public static final String EPOCH_SIZES = "populationSizes";
    public static final String TREE_MODEL = "populationTree";

    public String getParserName() {
        return PIECEWISE_POPULATION;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        Units.Type units = XMLUnits.Utils.getUnitsAttr(xo);

        XMLObject cxo = xo.getChild(EPOCH_SIZES);
        Parameter epochSizes = (Parameter) cxo.getChild(Parameter.class);

        cxo = xo.getChild(TREE_MODEL);
        TreeModel treeModel = (TreeModel) cxo.getChild(TreeModel.class);

        boolean isLinear = xo.getBooleanAttribute("linear");

        return new ScaledPiecewiseModel(epochSizes, treeModel, isLinear, units);
    }


    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "This element represents a piecewise population model";
    }

    public Class getReturnType() {
        return PiecewisePopulationModel.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(EPOCH_SIZES, new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
            new ElementRule(TREE_MODEL, new XMLSyntaxRule[]{new ElementRule(TreeModel.class)}),
            XMLUnits.SYNTAX_RULES[0],
            AttributeRule.newBooleanRule("linear")
    };

}
