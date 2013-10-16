package dr.evomodelxml.coalescent;

import dr.evomodel.coalescent.PopulationSizeGraph;
import dr.evomodel.coalescent.VariableDemographicModel;
import dr.xml.*;

/**
 */
public class PopulationSizeGraphParser extends AbstractXMLObjectParser {

    public static String POPGRAPH_STATISTIC = "popGraph";

    public String getParserName() {
        return POPGRAPH_STATISTIC;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        Object child = xo.getChild(0);

        if (child instanceof VariableDemographicModel) {
            final double dim = xo.getDoubleAttribute("time");
            return new PopulationSizeGraph((VariableDemographicModel) child, dim);
        }

        throw new XMLParseException("Unknown element found in " + getParserName() + " element:" + child);
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "This element returns a statistic that is the population size at evenly spaced intervals over tree.";
    }

    public Class getReturnType() {
        return PopulationSizeGraph.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
            new ElementRule(VariableDemographicModel.class, 1, 1),
            AttributeRule.newDoubleRule("time")
    };
}
