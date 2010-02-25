package dr.inferencexml.model;

import dr.inference.model.Statistic;
import dr.inference.model.SubStatistic;
import dr.xml.*;

/**
 */
public class SubStatisticParser extends AbstractXMLObjectParser {

    public static final String SUB_STATISTIC = "subStatistic";
    public static final String DIMENSION = "dimension";

    public String getParserName() {
        return SUB_STATISTIC;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        String name;
        if (xo.hasAttribute(Statistic.NAME) || xo.hasAttribute(dr.xml.XMLParser.ID))
            name = xo.getAttribute(Statistic.NAME, xo.getId());
        else
            name = "";

        final Statistic stat = (Statistic) xo.getChild(Statistic.class);

        final int[] values = xo.getIntegerArrayAttribute(DIMENSION);

        if (values.length == 0) {
            throw new XMLParseException("Must specify at least one dimension");
        }

        final int dim = stat.getDimension();

        for (int value : values) {
            if (value >= dim || value < 0) {
                throw new XMLParseException("Dimension " + value + " is not a valid dimension.");
            }
        }

        return new SubStatistic(name, values, stat);
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "Allows you to choose specific dimensions of a given statistic";
    }

    public Class getReturnType() {
        return SubStatistic.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return new XMLSyntaxRule[]{
                AttributeRule.newIntegerArrayRule(DIMENSION, false),
                new ElementRule(Statistic.class),
        };
    }
}
