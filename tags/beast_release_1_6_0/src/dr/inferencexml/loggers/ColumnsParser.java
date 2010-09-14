package dr.inferencexml.loggers;

import dr.inference.loggers.*;
import dr.util.Identifiable;
import dr.xml.*;

import java.util.ArrayList;

/**
 *
 */
public class ColumnsParser extends AbstractXMLObjectParser {
    public static final String COLUMN = "column";
    public static final String LABEL = "label";
    public static final String SIGNIFICANT_FIGURES = "sf";
    public static final String DECIMAL_PLACES = "dp";
    public static final String WIDTH = "width";
    public static final String FORMAT = "format";
    public static final String PERCENT = "percent";
    public static final String BOOL = "boolean";

    public String getParserName() {
        return COLUMN;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        String label = xo.getAttribute(LABEL, "");
        final int sf = xo.getAttribute(SIGNIFICANT_FIGURES, -1);
        final int dp = xo.getAttribute(DECIMAL_PLACES, -1);
        final int width = xo.getAttribute(WIDTH, -1);

        String format = xo.getAttribute(FORMAT, "");

        ArrayList colList = new ArrayList();

        for (int i = 0; i < xo.getChildCount(); i++) {

            Object child = xo.getChild(i);
            LogColumn[] cols;

            if (child instanceof Loggable) {
                cols = ((Loggable) child).getColumns();
            } else if (child instanceof Identifiable) {
                cols = new LogColumn[]{new LogColumn.Default(((Identifiable) child).getId(), child)};
            } else {
                cols = new LogColumn[]{new LogColumn.Default(child.getClass().toString(), child)};
            }

            if (format.equals(PERCENT)) {
                for (int k = 0; k < cols.length; ++k) {
                    if (cols[k] instanceof NumberColumn) {
                        cols[k] = new PercentColumn((NumberColumn) cols[k]);
                    }
                }
            } else if (format.equals(BOOL)) {
                for (int k = 0; k < cols.length; ++k) {
                    if (cols[k] instanceof NumberColumn) {
                        cols[k] = new BooleanColumn((NumberColumn) cols[k]);
                    }
                }
            }

            for (int j = 0; j < cols.length; j++) {

                if (!label.equals("")) {
                    if (cols.length > 1) {
                        cols[j].setLabel(label + Integer.toString(j + 1));
                    } else {
                        cols[j].setLabel(label);
                    }
                }

                if (cols[j] instanceof NumberColumn) {
                    if (sf != -1) {
                        ((NumberColumn) cols[j]).setSignificantFigures(sf);
                    }
                    if (dp != -1) {
                        ((NumberColumn) cols[j]).setDecimalPlaces(dp);
                    }
                }

                if (width > 0) {
                    cols[j].setMinimumWidth(width);
                }

                colList.add(cols[j]);
            }
        }

        LogColumn[] columns = new LogColumn[colList.size()];
        colList.toArray(columns);

        return new Columns(columns);
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "Specifies formating options for one or more columns in a log file.";
    }

    public Class getReturnType() {
        return Columns.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
            new StringAttributeRule(LABEL,
                    "The label of the column. " +
                            "If this is specified and more than one statistic is in this column, " +
                            "then the label will be appended by the index of the statistic to create individual column names", true),
            AttributeRule.newIntegerRule(SIGNIFICANT_FIGURES, true),
            AttributeRule.newIntegerRule(DECIMAL_PLACES, true),
            AttributeRule.newIntegerRule(WIDTH, true),
            AttributeRule.newStringRule(FORMAT, true),
            // Anything goes???
            new ElementRule(Object.class, 1, Integer.MAX_VALUE),
    };

}
