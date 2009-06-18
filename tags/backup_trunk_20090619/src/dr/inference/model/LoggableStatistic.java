package dr.inference.model;

import dr.xml.*;
import dr.inference.loggers.Loggable;
import dr.inference.loggers.LogColumn;
import dr.inference.loggers.NumberColumn;

/**
 * @author Marc A. Suchard
 */
public class LoggableStatistic extends Statistic.Abstract {

    public static final String LOGGABLE_STATISTIC = "loggableStatistic";

    private final NumberColumn[] logColumns;

    public LoggableStatistic(String name, NumberColumn[] columns) {
        super(name);
        logColumns = columns;
    }

    public int getDimension() {
        return logColumns.length;
    }

    public double getStatisticValue(int dim) {
        return logColumns[dim].getDoubleValue();
    }

    public String getDimensionName(int dim) {
        return logColumns[dim].getLabel();
    }

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return LOGGABLE_STATISTIC;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            String name;
            if (xo.hasAttribute(NAME) || xo.hasAttribute(ID))
                name = xo.getAttribute(NAME, xo.getId());
            else
                name = "";
                        

            final Loggable loggable = (Loggable) xo.getChild(Loggable.class);
            final LogColumn[] logColumns = loggable.getColumns();
            final NumberColumn[] numberColumns = new NumberColumn[logColumns.length];

            for(int i=0; i<logColumns.length; i++) {
                if (logColumns[i] instanceof NumberColumn) {
                    numberColumns[i] = (NumberColumn) logColumns[i];
                } else
                    throw new XMLParseException("Can only convert NumberColumns into Statistics");
            }

            return new LoggableStatistic(name, numberColumns);
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
                    new ElementRule(Loggable.class),
            };
        }
    };
}
