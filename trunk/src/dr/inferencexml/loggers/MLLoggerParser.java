package dr.inferencexml.loggers;

import dr.inference.loggers.*;
import dr.inference.model.Likelihood;
import dr.util.Identifiable;
import dr.xml.*;

import java.io.PrintWriter;

/**
 *
 */
public class MLLoggerParser extends LoggerParser {

    public static final String LOG_ML = "logML";
    public static final String LIKELIHOOD = "ml";

    public String getParserName() {
        return LOG_ML;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        Likelihood likelihood = (Likelihood) xo.getElementFirstChild(LIKELIHOOD);

        // logEvery of zero only displays at the end
        int logEvery = xo.getAttribute(LOG_EVERY, 0);

        final PrintWriter pw = getLogFile(xo, getParserName());

        LogFormatter formatter = new TabDelimitedFormatter(pw);

        MLLogger logger = new MLLogger(likelihood, formatter, logEvery);

        if (xo.hasAttribute(TITLE)) {
            logger.setTitle(xo.getStringAttribute(TITLE));
        }

        for (int i = 0; i < xo.getChildCount(); i++) {
            Object child = xo.getChild(i);

            if (child instanceof Columns) {

                logger.addColumns(((Columns) child).getColumns());

            } else if (child instanceof Loggable) {

                logger.add((Loggable) child);

            } else if (child instanceof Identifiable) {

                logger.addColumn(new LogColumn.Default(((Identifiable) child).getId(), child));

            } else {

                logger.addColumn(new LogColumn.Default(child.getClass().toString(), child));
            }
        }

        return logger;
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            AttributeRule.newIntegerRule(LOG_EVERY, true),
            new ElementRule(LIKELIHOOD,
                    new XMLSyntaxRule[]{new ElementRule(Likelihood.class)}),
            new OrRule(
                    new ElementRule(Columns.class, 1, Integer.MAX_VALUE),
                    new ElementRule(Loggable.class, 1, Integer.MAX_VALUE)
            )
    };

    public String getParserDescription() {
        return "Logs one or more items every time the given likelihood improves";
    }

    public Class getReturnType() {
        return MLLogger.class;
    }
}
