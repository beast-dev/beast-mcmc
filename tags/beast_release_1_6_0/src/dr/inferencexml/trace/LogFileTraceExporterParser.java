package dr.inferencexml.trace;

import dr.inference.trace.LogFileTraceExporter;
import dr.util.FileHelpers;
import dr.xml.*;

import java.io.File;

/**
 *
 */
public class LogFileTraceExporterParser extends AbstractXMLObjectParser {
    public static final String LOG_FILE_TRACE = "logFileTrace";

    private static final String FILENAME = "fileName";
    private static final String BURN_IN = "burnIn";

    public String getParserName() {
        return LOG_FILE_TRACE;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        final File file = FileHelpers.getFile(xo.getStringAttribute(FILENAME));
        final int burnIn = xo.getAttribute(BURN_IN, -1);

        try {
            return new LogFileTraceExporter(file, burnIn);
        } catch (Exception e) {
            throw new XMLParseException(e.getMessage());
        }
    }

    public String getParserDescription() {
        return "reconstruct population graph from variable dimension run.";
    }

    public Class getReturnType() {
        return LogFileTraceExporter.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            AttributeRule.newStringRule(FILENAME, false, "trace log."),
            AttributeRule.newIntegerRule(BURN_IN, true,
                    "The number of states (not sampled states, but actual states) that are discarded from the" +
                            " beginning of the trace before doing the analysis"),
    };
}
