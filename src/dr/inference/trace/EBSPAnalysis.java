package dr.inference.trace;

import dr.evolution.io.Importer;
import dr.evolution.io.NexusImporter;
import dr.evolution.io.TreeImporter;
import dr.evolution.tree.Tree;
import dr.evomodel.coalescent.VDdemographicFunction;
import dr.evomodel.coalescent.VariableDemographicModel;
import dr.inference.loggers.MCLogger;
import dr.stats.DiscreteStatistics;
import dr.util.HeapSort;
import dr.util.TabularData;
import dr.xml.*;

import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.Arrays;

/**
 @author Joseph Heled
 */
public class EBSPAnalysis extends TabularData {

    private double[] xPoints;
    private double[] means;
    private double[] medians ;
    private double[][] hpdLower;
    private double[][] hpdHigh;
    private double[] HPDLevels;
    // each bin covers xPoints[-1]/coalBins.length
    private int[] coalBins;

    private boolean quantiles;

    EBSPAnalysis(File log, File[] treeFiles, VariableDemographicModel.Type modelType,
                                String firstColumnName, String firstIndicatorColumnName,
                                String rootHeightColumnName, int coalPointBins,
                                double burnIn,
                                double[] inHPDLevels, boolean quantiles, boolean logSpace, boolean  mid)
            throws IOException, Importer.ImportException, TraceException {

        LogFileTraces ltraces = new LogFileTraces(log.getCanonicalPath(), log);
        ltraces.loadTraces();
        ltraces.setBurnIn(0);
        final int runLengthIncludingBurnin = ltraces.getStateCount();

        int intBurnIn =  (int)Math.floor(burnIn < 1 ? runLengthIncludingBurnin * burnIn : burnIn);
        final int nStates = runLengthIncludingBurnin - intBurnIn;
        //intBurnIn *= ltraces.getStepSize();
        ltraces.setBurnIn(intBurnIn * ltraces.getStepSize());

        assert ltraces.getStateCount() == nStates;

        this.quantiles = quantiles;
        HPDLevels = (inHPDLevels != null ) ? inHPDLevels :  new double[]{0.95};

        int populationFirstColumn = -1;
        int indicatorsFirstColumn = -1;
        int rootHeightColumn = -1;

        for(int n = 0; n < ltraces.getTraceCount(); ++n) {
            final String traceName = ltraces.getTraceName(n);
            if( traceName.equals(firstColumnName) ) {
                populationFirstColumn = n;
            } else if( traceName.equals(firstIndicatorColumnName) ) {
                indicatorsFirstColumn = n;
            } else if( rootHeightColumnName != null && traceName.equals(rootHeightColumnName) ) {
                rootHeightColumn = n;
            }
        }

        if( populationFirstColumn < 0 || indicatorsFirstColumn < 0 ) {
            throw new TraceException("incorrect trace column names: unable to find populations/indicators");
        }

        double binSize = 0;
        if( coalPointBins > 0 ) {
            if( rootHeightColumn < 0 ) {
                throw new TraceException("incorrect tree height column");
            }
            double hSum = -0;
            double[] h = new double[1];
            for(int ns = 0; ns < nStates; ++ns) {
                ltraces.getStateValues(ns, h, rootHeightColumn);
                hSum += h[0];
            }
            binSize = hSum / (nStates * coalPointBins);
            coalBins = new int[coalPointBins];
            Arrays.fill(coalBins, 0);
        }

        TreeImporter[] treeImporters = new TreeImporter[treeFiles.length];
        final boolean isStepWise = modelType == VariableDemographicModel.Type.STEPWISE;

        int nIndicators = 0;

        for(int k = 0; k < treeFiles.length; ++k) {

           // System.err.println("burnin " + treeFiles[k] + "(" + k + ")");

            treeImporters[k] = new NexusImporter(new FileReader(treeFiles[k]));
            assert intBurnIn > 0;
            for(int z = 0; z < intBurnIn-1; ++z) {
                treeImporters[k].importNextTree();
            }
            nIndicators += treeImporters[k].importNextTree().getExternalNodeCount() - 1;
        }

        if( isStepWise ) {
            nIndicators -= 1;
        }

        final int nXaxisPoints = nIndicators + (isStepWise ? 1 : 0) + 1;
        xPoints = new double[nXaxisPoints];
        Arrays.fill(xPoints, 0.0);

        VDdemographicFunction[] allDemog = new VDdemographicFunction[nStates];
        {
            double[] indicators = new double[nIndicators];
            double[] pop = new double[nIndicators+1];
            Tree[] tt = new Tree[treeFiles.length];

            for(int ns = 0; ns < nStates; ++ns) {
                ltraces.getStateValues(ns, indicators, indicatorsFirstColumn);
                ltraces.getStateValues(ns, pop, populationFirstColumn);

                for(int nt = 0; nt < tt.length; ++nt) {
                    tt[nt] = treeImporters[nt].importNextTree();
                }
                final VDdemographicFunction demoFunction =
                        new VDdemographicFunction(tt, modelType, indicators, pop, logSpace, mid);
                double[] xs = demoFunction.allTimePoints();
                for(int k = 0; k < xs.length; ++k) {
                    //assert xs[k] >= 0 : " " + k + " " + xs[k];
                    xPoints[k+1] += xs[k];
                }
                if( coalPointBins > 0 ) {
                    for (double x : xs) {
                        coalBins[Math.min((int) (x / binSize), coalBins.length-1)]++;
                    }
                }
                allDemog[ns] = demoFunction;

                demoFunction.freeze();
            }
            for(int k = 0; k < xPoints.length; ++k) {
                xPoints[k] /= nStates;
            }
        }

        double[] popValues = new double[nStates];
        means = new double[nXaxisPoints];
        medians = new double[nXaxisPoints];
        hpdLower = new double[HPDLevels.length][];
        hpdHigh = new double[HPDLevels.length][];

        for(int i = 0; i < HPDLevels.length; ++i) {
            hpdLower[i] = new double[nXaxisPoints];
            hpdHigh[i] = new double[nXaxisPoints];
        }

        for(int nx = 0; nx < xPoints.length; ++nx) {
            final double x = xPoints[nx];

            for(int ns = 0; ns < nStates; ++ns) {
                popValues[ns] = allDemog[ns].getDemographic(x);
            }
            int[] indices = new int[popValues.length];
            HeapSort.sort(popValues, indices);

            means[nx] = DiscreteStatistics.mean(popValues);
            for(int i = 0; i < HPDLevels.length; ++i) {

                if( quantiles ) {
                    hpdLower[i][nx] = DiscreteStatistics.quantile((1-HPDLevels[i])/2, popValues, indices);
                    hpdHigh[i][nx] = DiscreteStatistics.quantile((1+HPDLevels[i])/2, popValues, indices);
                } else {
                    final double[] hpd = DiscreteStatistics.HPDInterval(HPDLevels[i], popValues, indices);
                    hpdLower[i][nx] = hpd[0];
                    hpdHigh[i][nx] = hpd[1];
                }
            }
            medians[nx] = DiscreteStatistics.median(popValues, indices);
//            for(int k = 0; k < indices.length; ++k) {
//                System.out.print(popValues[indices[k]]);  System.out.print(",");
//            }
        }
    }

    private final String[] columnNames = {"time", "mean", "median"};

    public int nColumns() {
        return columnNames.length + 2*HPDLevels.length + (coalBins != null ? 1 : 0);
    }

    public String columnName(int nColumn) {
        final int fixed = columnNames.length;
        if( nColumn < fixed) {
            return columnNames[nColumn];
        }
        nColumn -= fixed;
        if( nColumn < 2*HPDLevels.length) {
          final double p = HPDLevels[nColumn/2];
          final String s = ( nColumn % 2 == 0 ) ? "lower" : "upper";
          return (quantiles ? "cpd " : "hpd ") + s + " " + Math.round(p*100);
        }
        assert (nColumn - 2*HPDLevels.length) == 0;
        return "bins";
    }

    public int nRows() {
        return Math.max(xPoints.length, (coalBins != null ? coalBins.length : 0));
    }

    public Object data(int nRow, int nColumn) {
        switch( nColumn ) {
            case 0:
            {
                if( nRow < xPoints.length ) {
                    return xPoints[nRow];
                }
                break;
            }
            case 1:
            {
                if( nRow < means.length ) {
                    return means[nRow];
                }
                break;
            }
            case 2:
            {  if( nRow < medians.length ) {
                return medians[nRow];
            }
                break;
            }
            default:
            {
                final int j = nColumn - columnNames.length;
                if( j < 2*HPDLevels.length ) {
                    if( nRow < xPoints.length ) {
                        final int k = j /2;
                        if( 0 <= k && k < HPDLevels.length ) {
                            if( j % 2 == 0 ) {
                                return hpdLower[k][nRow];
                            } else {
                                return hpdHigh[k][nRow];
                            }
                        }
                    }
                } else {
                    if( nRow < coalBins.length ) {
                        return coalBins[nRow];
                    }
                }
                break;
            }
        }
        return "";
    }

    // should be local to PARSER, but this makes them non-accesible from the outside in Java.

    public static final String VD_ANALYSIS = "VDAnalysis";
    public static final String FILE_NAME = "fileName";
    public static final String BURN_IN = "burnIn";
    public static final String HPD_LEVELS = "Confidencelevels";
    public static final String QUANTILES = "useQuantiles";
    public static final String LOG_SPACE = VariableDemographicModel.LOG_SPACE;
    public static final String USE_MIDDLE = VariableDemographicModel.USE_MIDPOINTS;

    public static final String TREE_LOG = "treeOfLoci";

    public static final String LOG_FILE_NAME = "logFileName";
    public static final String TREE_FILE_NAMES = "treeFileNames";
    public static final String MODEL_TYPE = "populationModelType";
    public static final String POPULATION_FIRST_COLUMN = "populationFirstColumn";
    public static final String INDICATORS_FIRST_COLUMN = "indicatorsFirstColumn";
    public static final String ROOTHEIGHT_COLUMN = "rootheightColumn";
    public static final String NBINS = "nBins";

    public static dr.xml.XMLObjectParser PARSER = new dr.xml.AbstractXMLObjectParser() {

        private String getElementText(XMLObject xo, String childName) throws XMLParseException {
            return ((XMLObject)xo.getChild(childName)).getStringChild(0);
        }

        public String getParserName() {
            return VD_ANALYSIS;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            try {

                // 10% is brun-in default
                double burnin = 0.1;
                if (xo.hasAttribute(BURN_IN)) {
                    burnin = xo.getDoubleAttribute(BURN_IN);
                    if( burnin < 0 )  throw new XMLParseException("burnIn should be either between 0 and 1 or a positive number");
                }

                final double[] hpdLevels = xo.hasAttribute(HPD_LEVELS) ? xo.getDoubleArrayAttribute(HPD_LEVELS) : null;

                final File log = MCLogger.getFile( getElementText(xo, LOG_FILE_NAME) );
                final XMLObject treeFileNames = (XMLObject) xo.getChild(TREE_FILE_NAMES);
                final int nTrees = treeFileNames.getChildCount();
                File[] treeFiles = new File[nTrees];
                for(int k = 0; k < nTrees; ++k) {
                    treeFiles[k] = MCLogger.getFile( ((XMLObject)treeFileNames.getChild(k)).getStringChild(0) );
                }

                String modelTypeName = getElementText(xo, MODEL_TYPE).trim().toUpperCase();
                String populationFirstColumn = getElementText(xo, POPULATION_FIRST_COLUMN);
                String indicatorsFirstColumn = getElementText(xo, INDICATORS_FIRST_COLUMN);
                VariableDemographicModel.Type modelType = VariableDemographicModel.Type.valueOf(modelTypeName);

                String rootHeightColumn = null;
                int nBins = -1;
                if( xo.hasAttribute(NBINS) ) {
                    if( xo.getChild(ROOTHEIGHT_COLUMN) != null )  {
                        rootHeightColumn = getElementText(xo,ROOTHEIGHT_COLUMN);
                        nBins = xo.getIntegerAttribute(NBINS);
                    }
                }

                final boolean quantiles = xo.hasAttribute(QUANTILES) && xo.getBooleanAttribute(QUANTILES);
                final boolean logSpace = xo.hasAttribute(LOG_SPACE) && xo.getBooleanAttribute(LOG_SPACE);
                final boolean useMid = xo.hasAttribute(USE_MIDDLE) && xo.getBooleanAttribute(USE_MIDDLE);

                return new EBSPAnalysis(log, treeFiles, modelType,
                        populationFirstColumn, indicatorsFirstColumn,
                        rootHeightColumn, nBins,
                        burnin, hpdLevels, quantiles, logSpace, useMid);

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

        private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
                AttributeRule.newDoubleRule(BURN_IN, true, "The number of states (not sampled states, but" +
                        " actual states) that are discarded from the beginning of the trace and are excluded from " +
                        "the analysis"),
                AttributeRule.newDoubleArrayRule(HPD_LEVELS, true),
                AttributeRule.newBooleanRule(QUANTILES, true),
                AttributeRule.newBooleanRule(LOG_SPACE, true),
                AttributeRule.newBooleanRule(USE_MIDDLE, true),
                AttributeRule.newIntegerRule(NBINS, true),
                
                new ElementRule(LOG_FILE_NAME, String.class, "The name of a BEAST log file"),
                new ElementRule(TREE_FILE_NAMES,
                        new XMLSyntaxRule[] {
                             new ElementRule(TREE_LOG, String.class, "The name of a BEAST trees log file", 1, Integer.MAX_VALUE)
                        } ),
                new ElementRule(MODEL_TYPE, String.class, "population model type (stepwise, linear ..."),
                new ElementRule(POPULATION_FIRST_COLUMN, String.class, "Name of first column of population size"),
                new ElementRule(INDICATORS_FIRST_COLUMN, String.class, "Name of first column of population indicators"),
                new ElementRule(ROOTHEIGHT_COLUMN, String.class, "Name of trace column of root height", true),
        };
    };
}
