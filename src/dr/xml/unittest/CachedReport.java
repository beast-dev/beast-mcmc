package dr.xml.unittest;

import dr.xml.*;

public class CachedReport implements Reportable {
    private String report;
    private final Reportable reportable;

    CachedReport(Reportable reportable) {
        this.reportable = reportable;
        this.report = null;
    }

    @Override
    public String getReport() {
        if (report == null) {
            report = reportable.getReport();
        }
        return report;
    }

    private static final String CACHED_REPORT = "cachedReport";


    public static final AbstractXMLObjectParser PARSER = new AbstractXMLObjectParser() {
        @Override
        public Class getReturnType() {
            return CachedReport.class;
        }

        @Override
        public Object parseXMLObject(XMLObject xo) throws XMLParseException {
            Reportable reportable = (Reportable) xo.getChild(Reportable.class);
            return new CachedReport(reportable);
        }

        @Override
        public XMLSyntaxRule[] getSyntaxRules() {
            return new XMLSyntaxRule[]{
                    new ElementRule(Reportable.class)
            };
        }

        @Override
        public String getParserDescription() {
            return "Reportable object that caches the report rather than re-computing it.";
        }

        @Override
        public String getParserName() {
            return CACHED_REPORT;
        }


    };
}
