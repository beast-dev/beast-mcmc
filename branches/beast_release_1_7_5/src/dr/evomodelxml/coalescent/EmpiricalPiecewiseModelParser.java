package dr.evomodelxml.coalescent;

import dr.evolution.util.Units;
import dr.evomodel.coalescent.EmpiricalPiecewiseModel;
import dr.evoxml.util.XMLUnits;
import dr.inference.model.Parameter;
import dr.xml.*;

/**
 * Parses an element from an DOM document into a PiecewisePopulation.
 */
public class EmpiricalPiecewiseModelParser extends AbstractXMLObjectParser {

    public static final String EMPIRICAL_PIECEWISE = "empiricalPiecewise";
    public static final String INTERVAL_WIDTHS = "intervalWidths";
    public static final String POPULATION_SIZES = "populationSizes";
    public static final String TAU = "generationLength";
    public static final String THRESHOLD = "threshold";
    public static final String LAG = "lag";

    public String getParserName() {
        return EMPIRICAL_PIECEWISE;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        Units.Type units = XMLUnits.Utils.getUnitsAttr(xo);


        XMLObject cxo = xo.getChild(INTERVAL_WIDTHS);
        double[] intervalWidths = cxo.getDoubleArrayAttribute("values");

        cxo = xo.getChild(POPULATION_SIZES);
        Parameter popSizes = (Parameter) cxo.getChild(Parameter.class);

        cxo = xo.getChild(TAU);
        Parameter scaleParam = (Parameter) cxo.getChild(Parameter.class);

        cxo = xo.getChild(THRESHOLD);
        Parameter bParam = (Parameter) cxo.getChild(Parameter.class);

        cxo = xo.getChild(LAG);
        Parameter offsetParam = (Parameter) cxo.getChild(Parameter.class);

        return new EmpiricalPiecewiseModel(intervalWidths, popSizes, scaleParam, bParam, offsetParam, units);
    }


    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "This element represents a piecewise population model";
    }

    public Class getReturnType() {
        return EmpiricalPiecewiseModel.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(INTERVAL_WIDTHS,
                    new XMLSyntaxRule[]{AttributeRule.newDoubleArrayRule("values", false),}),
            XMLUnits.SYNTAX_RULES[0],
            new ElementRule(POPULATION_SIZES,
                    new XMLSyntaxRule[]{new ElementRule(Parameter.class)},
                    "The effective population sizes of each interval."),
            new ElementRule(TAU,
                    new XMLSyntaxRule[]{new ElementRule(Parameter.class)},
                    "The scale factor."),
            new ElementRule(THRESHOLD,
                    new XMLSyntaxRule[]{new ElementRule(Parameter.class)},
                    "The threshold before counts occur."),
            new ElementRule(LAG,
                    new XMLSyntaxRule[]{new ElementRule(Parameter.class)},
                    "The lag between actual population sizes and genetic diversity.")
    };
}
