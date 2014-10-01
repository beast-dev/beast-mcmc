package dr.evomodelxml.coalescent;

import dr.evolution.util.Units;
import dr.evomodel.coalescent.PiecewisePopulationModel;
import dr.evoxml.util.XMLUnits;
import dr.inference.model.Parameter;
import dr.xml.*;

/**
 * Parses an element from an DOM document into a PiecewisePopulation.
 */
public class PiecewisePopulationModelParser extends AbstractXMLObjectParser {

    public static final String PIECEWISE_POPULATION = "piecewisePopulation";
    public static final String EPOCH_SIZES = "epochSizes";
    public static final String POPULATION_SIZE = "populationSize";
    public static final String GROWTH_RATES = "growthRates";
    public static final String EPOCH_WIDTHS = "epochWidths";

    public static final String WIDTHS = "widths";
    public static final String LINEAR = "linear";

    public String getParserName() {
        return PIECEWISE_POPULATION;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        Units.Type units = XMLUnits.Utils.getUnitsAttr(xo);

        XMLObject obj = xo.getChild(EPOCH_WIDTHS);
        double[] epochWidths = obj.getDoubleArrayAttribute(WIDTHS);

        if (xo.hasChildNamed(EPOCH_SIZES)) {
            Parameter epochSizes = (Parameter) xo.getElementFirstChild(EPOCH_SIZES);

            boolean isLinear = false;
            if (xo.hasAttribute(LINEAR)) {
                isLinear = xo.getBooleanAttribute(LINEAR);
            }

            return new PiecewisePopulationModel(PIECEWISE_POPULATION, epochSizes, epochWidths, isLinear, units);
        } else {
            Parameter populationSize = (Parameter) xo.getElementFirstChild(POPULATION_SIZE);
            Parameter growthRates = (Parameter) xo.getElementFirstChild(GROWTH_RATES);
            return new PiecewisePopulationModel(PIECEWISE_POPULATION, populationSize, growthRates, epochWidths, units);
        }
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

    private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
            new XORRule(
                    new ElementRule(EPOCH_SIZES,
                            new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
                    new AndRule(
                            new ElementRule(POPULATION_SIZE,
                                    new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
                            new ElementRule(GROWTH_RATES,
                                    new XMLSyntaxRule[]{new ElementRule(Parameter.class)})
                    )
            ),
            new ElementRule(EPOCH_WIDTHS,
                    new XMLSyntaxRule[]{AttributeRule.newDoubleArrayRule(WIDTHS)}),
            AttributeRule.newBooleanRule(LINEAR, true)
    };
}
