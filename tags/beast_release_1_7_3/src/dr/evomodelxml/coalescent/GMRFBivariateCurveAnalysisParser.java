package dr.evomodelxml.coalescent;

import dr.evomodel.coalescent.GMRFBivariateCurveAnalysis;
import dr.inferencexml.trace.TraceAnalysisParser;
import dr.xml.*;

import java.io.File;

/**
 *
 */
public class GMRFBivariateCurveAnalysisParser extends AbstractXMLObjectParser {

    private static final String GMRF_BIVARIATE_CURVE_ANALYSIS = "gmrfBivariateCurveAnalysis";
    private static final String END_TIME = "endTime";
    private static final String FILE_NAME_ONE = "fileName1";
    private static final String FILE_NAME_TWO = "fileName2";

    public String getParserName() {
        return GMRF_BIVARIATE_CURVE_ANALYSIS;
    }

    public String getParserDescription() {
        return "Integrates two curves";
    }

    public Class getReturnType() {
        return GMRFBivariateCurveAnalysis.class;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        String[] inputFileName = {
                xo.getStringAttribute(FILE_NAME_ONE),
                xo.getStringAttribute(FILE_NAME_TWO)};

        File[] file = new File[2];
        String[] name = new String[2];
        String[] parent = new String[2];

        for (int i = 0; i < inputFileName.length; i++) {
            file[i] = new File(inputFileName[i]);
            name[i] = file[i].getName();
            parent[i] = file[i].getParent();

            if (!file[i].isAbsolute()) {
                parent[i] = System.getProperty("user.dir");
            }

            file[i] = new File(parent[i], name[i]);
            inputFileName[i] = file[i].getName();
        }

        int burnin = -1;
        if (xo.hasAttribute(TraceAnalysisParser.BURN_IN)) {
            burnin = xo.getIntegerAttribute(TraceAnalysisParser.BURN_IN);
        }

        double endTime = xo.getDoubleAttribute(END_TIME);

        try {
            GMRFBivariateCurveAnalysis analysis = new GMRFBivariateCurveAnalysis(inputFileName, endTime, burnin);

//            	analysis.integrationReport();
            analysis.densityReport();
        } catch (Exception e) {
            System.err.println(e.getMessage());
        }
        return null;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return new XMLSyntaxRule[]{
                AttributeRule.newStringRule(FILE_NAME_ONE, false),
                AttributeRule.newStringRule(FILE_NAME_TWO, false),
                AttributeRule.newStringRule(END_TIME, false),
                AttributeRule.newIntegerRule(TraceAnalysisParser.BURN_IN, true),
        };
    }
    
}
