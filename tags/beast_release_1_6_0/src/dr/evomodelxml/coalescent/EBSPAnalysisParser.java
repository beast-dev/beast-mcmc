package dr.evomodelxml.coalescent;

import dr.evolution.io.Importer;
import dr.evomodel.coalescent.EBSPAnalysis;
import dr.evomodel.coalescent.VariableDemographicModel;
import dr.inference.trace.TraceException;
import dr.util.FileHelpers;
import dr.xml.*;

import java.io.File;
import java.io.FileWriter;
import java.io.PrintWriter;

/**
 */
public class EBSPAnalysisParser extends AbstractXMLObjectParser {

    public static final String VD_ANALYSIS = "VDAnalysis";
    public static final String FILE_NAME = "fileName";
    public static final String BURN_IN = "burnIn";
    public static final String HPD_LEVELS = "Confidencelevels";
    public static final String QUANTILES = "useQuantiles";
    public static final String LOG_SPACE = VariableDemographicModelParser.LOG_SPACE;
    public static final String USE_MIDDLE = VariableDemographicModelParser.USE_MIDPOINTS;
    public static final String N_CHANGES = "nChanges";

    public static final String TREE_LOG = "treeOfLoci";

    public static final String LOG_FILE_NAME = "logFileName";
    public static final String TREE_FILE_NAMES = "treeFileNames";
    public static final String MODEL_TYPE = "populationModelType";
    public static final String POPULATION_FIRST_COLUMN = "populationFirstColumn";
    public static final String INDICATORS_FIRST_COLUMN = "indicatorsFirstColumn";
    public static final String ROOTHEIGHT_COLUMN = "rootheightColumn";
    public static final String ALLDEMO_COLUMN = "allDemographicsFileName";
    public static final String NBINS = "nBins";

    private String getElementText(XMLObject xo, String childName) throws XMLParseException {
        return xo.getChild(childName).getStringChild(0);
    }

    public String getParserName() {
        return VD_ANALYSIS;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        try {

            // 10% is brun-in default
            final double burnin = xo.getAttribute(BURN_IN, 0.1);
            if (burnin < 0) {
                throw new XMLParseException("burnIn should be either between 0 and 1 or a positive number");
            }

            final double[] hpdLevels = xo.hasAttribute(HPD_LEVELS) ? xo.getDoubleArrayAttribute(HPD_LEVELS) : null;

            final File log = FileHelpers.getFile(getElementText(xo, LOG_FILE_NAME));
            final XMLObject treeFileNames = xo.getChild(TREE_FILE_NAMES);
            final int nTrees = treeFileNames.getChildCount();
            File[] treeFiles = new File[nTrees];
            for (int k = 0; k < nTrees; ++k) {
                treeFiles[k] = FileHelpers.getFile(((XMLObject) treeFileNames.getChild(k)).getStringChild(0));
            }

            String modelTypeName = getElementText(xo, MODEL_TYPE).trim().toUpperCase();
            String populationFirstColumn = getElementText(xo, POPULATION_FIRST_COLUMN);
            String indicatorsFirstColumn = getElementText(xo, INDICATORS_FIRST_COLUMN);
            VariableDemographicModel.Type modelType = VariableDemographicModel.Type.valueOf(modelTypeName);

            String rootHeightColumn = null;
            int nBins = -1;
            if (xo.hasAttribute(NBINS)) {
                if (xo.getChild(ROOTHEIGHT_COLUMN) != null) {
                    rootHeightColumn = getElementText(xo, ROOTHEIGHT_COLUMN);
                    nBins = xo.getIntegerAttribute(NBINS);
                }
            }

            PrintWriter allDemoWriter = null;
            if (xo.getChild(ALLDEMO_COLUMN) != null) {
                String fName = getElementText(xo, ALLDEMO_COLUMN);
                allDemoWriter = new PrintWriter(new FileWriter(fName));
            }

            final boolean quantiles = xo.getAttribute(QUANTILES, false);
            final boolean logSpace = xo.getAttribute(LOG_SPACE, false);
            final boolean useMid = xo.getAttribute(USE_MIDDLE, false);
            final int onlyNchanges = xo.getAttribute(N_CHANGES, -1);

            return new EBSPAnalysis(log, treeFiles, modelType,
                    populationFirstColumn, indicatorsFirstColumn,
                    rootHeightColumn, nBins,
                    burnin, hpdLevels, quantiles, logSpace, useMid, onlyNchanges,
                    allDemoWriter);

        } catch (java.io.IOException ioe) {
            throw new XMLParseException(ioe.getMessage());
        } catch (Importer.ImportException e) {
            throw new XMLParseException(e.toString());
        } catch (TraceException e) {
            throw new XMLParseException(e.toString());
        }
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "reconstruct population graph from EBSP run.";
    }

    public Class getReturnType() {
        return EBSPAnalysis.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            AttributeRule.newDoubleRule(BURN_IN, true, "The number of states (not sampled states, but" +
                    " actual states) that are discarded from the beginning of the trace and are excluded from " +
                    "the analysis"),
            AttributeRule.newDoubleArrayRule(HPD_LEVELS, true),
            AttributeRule.newBooleanRule(QUANTILES, true),
            AttributeRule.newBooleanRule(LOG_SPACE, true),
            AttributeRule.newBooleanRule(USE_MIDDLE, true),
            AttributeRule.newIntegerRule(NBINS, true),
            AttributeRule.newIntegerRule(N_CHANGES, true),

            new ElementRule(LOG_FILE_NAME, String.class, "The name of a BEAST log file"),
            new ElementRule(TREE_FILE_NAMES,
                    new XMLSyntaxRule[]{
                            new ElementRule(TREE_LOG, String.class, "The name of a BEAST trees log file", 1, Integer.MAX_VALUE)
                    }),
            new ElementRule(MODEL_TYPE, String.class, "population model type (stepwise, linear ..."),
            new ElementRule(POPULATION_FIRST_COLUMN, String.class, "Name of first column of population size"),
            new ElementRule(INDICATORS_FIRST_COLUMN, String.class, "Name of first column of population indicators"),
            new ElementRule(ROOTHEIGHT_COLUMN, String.class, "Name of trace column of root height", true),
            new ElementRule(ALLDEMO_COLUMN, String.class, "Name of file to output all demographic functions", true),
    };
}
