package dr.inference.xml;

import dr.app.beast.BeastVersion;
import dr.inference.loggers.*;
import dr.math.MathUtils;
import dr.util.FileHelpers;
import dr.util.Identifiable;
import dr.util.Property;
import dr.xml.*;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.PrintWriter;
import java.util.Date;

/**
 * @author Alexei Drummond
 * @author Andrew Rambaut
 */
public class LoggerParser extends AbstractXMLObjectParser {

    public static final String LOG = "log";
    public static final String ECHO = "echo";
    public static final String ECHO_EVERY = "echoEvery";
    public static final String TITLE = "title";
    public static final String FILE_NAME = FileHelpers.FILE_NAME;
    public static final String FORMAT = "format";
    public static final String TAB = "tab";
    public static final String HTML = "html";
    public static final String PRETTY = "pretty";
    public static final String LOG_EVERY = "logEvery";

    public static final String COLUMNS = "columns";
    public static final String COLUMN = "column";
    public static final String LABEL = "label";
    public static final String SIGNIFICANT_FIGURES = "sf";
    public static final String DECIMAL_PLACES = "dp";
    public static final String WIDTH = "width";

    public String getParserName() {
        return LOG;
    }

    /**
     * @return an object based on the XML element it was passed.
     */
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        // You must say how often you want to log
        final int logEvery = xo.getIntegerAttribute(LOG_EVERY);

        String fileName = null;
        if (xo.hasAttribute(FILE_NAME)) {
            fileName = xo.getStringAttribute(FILE_NAME);
        }

        final PrintWriter pw = getLogFile(xo, getParserName());

        final LogFormatter formatter = new TabDelimitedFormatter(pw);

        final MCLogger logger = new MCLogger(fileName, formatter, logEvery, !xo.hasAttribute(FILE_NAME));

        if (xo.hasAttribute(TITLE)) {
            logger.setTitle(xo.getStringAttribute(TITLE));
        } else {

            final BeastVersion version = new BeastVersion();

            final String title = "BEAST " + version.getVersionString() +
                    ", " + version.getBuildString() + "\n" +

                    "Generated " + (new Date()).toString() + " [seed=" + MathUtils.getSeed() + "]";
            logger.setTitle(title);
        }

        for (int i = 0; i < xo.getChildCount(); i++) {

            final Object child = xo.getChild(i);

            if (child instanceof Columns) {

                logger.addColumns(((Columns) child).getColumns());

            } else if (child instanceof Loggable) {

                logger.add((Loggable) child);

            } else if (child instanceof Identifiable) {

                logger.addColumn(new LogColumn.Default(((Identifiable) child).getId(), child));

            } else if (child instanceof Property) {
                logger.addColumn(new LogColumn.Default(((Property) child).getAttributeName(), child));
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
            AttributeRule.newIntegerRule(LOG_EVERY),
            new StringAttributeRule(FILE_NAME,
                    "The name of the file to send log output to. " +
                            "If no file name is specified then log is sent to standard output", true),
            new StringAttributeRule(TITLE,
                    "The title of the log", true),
            new OrRule(
                    new XMLSyntaxRule[]{
                            new ElementRule(Columns.class, 1, Integer.MAX_VALUE),
                            new ElementRule(Loggable.class, 1, Integer.MAX_VALUE),
                            new ElementRule(Object.class, 1, Integer.MAX_VALUE)
                    }
            )
    };

    public String getParserDescription() {
        return "Logs one or more items at a given frequency to the screen or to a file";
    }

    public Class getReturnType() {
        return MLLogger.class;
    }

    /**
     * Allow a file relative to beast xml file with a prefix of ./
     *
     * @param xo         log element
     * @param parserName for error messages
     * @return logger object
     * @throws XMLParseException if file can't be created for some reason
     */
    public static PrintWriter getLogFile(XMLObject xo, String parserName) throws XMLParseException {
        if (xo.hasAttribute(FILE_NAME)) {

            final String fileName = xo.getStringAttribute(FILE_NAME);

            final File logFile = FileHelpers.getFile(fileName);

//	         System.out.println("Writing log file to "+parent+System.getProperty("path.separator")+name);
            try {
                return new PrintWriter(new FileOutputStream(logFile));
            } catch (FileNotFoundException fnfe) {
                throw new XMLParseException("File '" + logFile.getAbsolutePath() + "' can not be opened for " + parserName + " element.");
            }

        }
        return new PrintWriter(System.out);
    }
}

