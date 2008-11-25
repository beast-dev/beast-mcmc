package dr.evomodelxml;

import dr.app.beast.BeastVersion;
import dr.inference.loggers.*;
import dr.math.MathUtils;
import dr.util.Identifiable;
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
    public static final String FILE_NAME = "fileName";
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

    private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
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

            final File logFile = getFile(fileName);

//	         System.out.println("Writing log file to "+parent+System.getProperty("path.separator")+name);
            try {
                return new PrintWriter(new FileOutputStream(logFile));
            } catch (FileNotFoundException fnfe) {
                throw new XMLParseException("File '" + logFile.getAbsolutePath() + "' can not be opened for " + parserName + " element.");
            }

        }
        return new PrintWriter(System.out);
    }

    /**
     * Resolve file from name.
     * <p/>
     * Keep A fully qualified (i.e. absolute path) as is. A name starting with a "./" is
     * relative to the master BEAST directory. Any other name is stripped of any directory
     * component and placed in the "user.dir" directory.
     *
     * @param fileName an absolute or relative file name
     * @return a File object resolved from provided file name
     */
    public static File getFile(String fileName) {
        final boolean localFile = fileName.startsWith("./");
        final boolean relative = masterBeastDirectory != null && localFile;
        if (localFile) {
            fileName = fileName.substring(2);
        }

        final File file = new File(fileName);
        final String name = file.getName();
        String parent = file.getParent();

        if (!file.isAbsolute()) {
            String p;
            if (relative) {
                p = masterBeastDirectory.getAbsolutePath();
            } else {
                p = System.getProperty("user.dir");
            }
            if (parent != null && parent.length() > 0) {
                parent = p + '/' + parent;
            } else {
                parent = p;
            }
        }
        return new File(parent, name);
    }

    // directory where beast xml file resides
    private static File masterBeastDirectory = null;

    public static void setMasterDir(File fileName) {
        masterBeastDirectory = fileName;
    }
}

