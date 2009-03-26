package dr.inference.trace;

import dr.evomodelxml.LoggerParser;
import dr.util.TabularData;
import dr.xml.AttributeRule;
import dr.xml.XMLObject;
import dr.xml.XMLParseException;
import dr.xml.XMLSyntaxRule;

import java.io.File;
import java.io.IOException;

/**
 * Export trace analysis data such as mean,median,HPD and ESS of trace variables.
 *
 * @author joseph
 *         Date: 25/10/2007
 */
public class LogFileTraceExporter extends TabularData {
    private final LogFileTraces analysis;
    private final String[] rows = {"mean", "median", "hpdLower", "hpdUpper", "ESS"};
    TraceDistribution[] distributions;

    public LogFileTraceExporter(File file, int burnin) throws TraceException, IOException {

        analysis = new LogFileTraces(file.getCanonicalPath(), file);
        analysis.loadTraces();
        if (burnin >= 0) {
            analysis.setBurnIn(burnin);
        }

        distributions = new TraceDistribution[nColumns()];
    }

    public int nColumns() {
        return analysis.getTraceCount();
    }

    public String columnName(int nColumn) {
        return analysis.getTraceName(nColumn);
    }

    public int nRows() {
        return rows.length; // analysis.getStateCount();
    }

    public Object data(int nRow, int nColumn) {
        if (distributions[nColumn] == null) {
            analysis.analyseTrace(nColumn);
            distributions[nColumn] = analysis.getDistributionStatistics(nColumn);
        }

        TraceDistribution distribution = distributions[nColumn];

        switch (nRow) {
            case 0: {
                return distribution.getMean();
            }
            case 1: {
                return distribution.getMedian();
            }
            case 2: {
                return distribution.getLowerHPD();
            }
            case 3: {
                return distribution.getUpperHPD();
            }
            case 4: {
                return distribution.getESS();
            }
        }

        return null;
        // return analysis.getStateValue(nColumn,  nRow);
    }

    public static dr.xml.XMLObjectParser PARSER = new dr.xml.AbstractXMLObjectParser() {
        private static final String FILENAME = "fileName";
        private static final String BURN_IN = "burnIn";

        public String getParserName() {
            return "logFileTrace";
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            final File file = LoggerParser.getFile(xo.getStringAttribute(FILENAME));
            int burnIn = xo.getAttribute(BURN_IN, -1);

            try {
                return new LogFileTraceExporter(file, burnIn);
            } catch (Exception e) {
                return new XMLParseException(e.getMessage());
            }
        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

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
    };
}
