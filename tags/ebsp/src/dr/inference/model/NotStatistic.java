package dr.inference.model;

import dr.xml.*;

/**
 * Wraps a statistic and returns (1.0 - value) where value is the wrapped
 * statistic.
 *
 *
 * @version $Id: NegativeStatistic.java,v 1.2 2005/05/24 20:26:00 rambaut Exp $
 *
 * @author Alexei Drummond
 * @author Andrew Rambaut
 */
public class NotStatistic extends Statistic.Abstract {

    public static String NOT_STATISTIC = "notStatistic";

    private Statistic statistic = null;

    public NotStatistic(String name, Statistic statistic) {
        super(name);
        this.statistic = statistic;
    }

    public int getDimension() {
        return statistic.getDimension();
    }

    /** @return mean of contained statistics */
    public double getStatisticValue(int dim) {

        return 1.0 - statistic.getStatisticValue(dim);
    }

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public String[] getParserNames() { return new String[] { getParserName(), "not" }; }
        public String getParserName() { return NotStatistic.NOT_STATISTIC; }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            NotStatistic notStatistic;

            Object child = xo.getChild(0);
            if (child instanceof Statistic) {
                notStatistic = new NotStatistic(NotStatistic.NOT_STATISTIC, (Statistic)child);
            } else {
                throw new XMLParseException("Unknown element found in " + getParserName() + " element:" + child);
            }

            return notStatistic;
        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public String getParserDescription() {
            return "This element returns a statistic that is the element-wise inverse of the child statistic.";
        }

        public Class getReturnType() { return NotStatistic.class; }

        public XMLSyntaxRule[] getSyntaxRules() { return rules; }

        private XMLSyntaxRule[] rules = new XMLSyntaxRule[] {
            new ElementRule(Statistic.class, 1, 1 )
        };
    };
}
