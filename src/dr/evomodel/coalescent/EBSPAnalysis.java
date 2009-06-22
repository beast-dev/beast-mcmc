package dr.evomodel.coalescent;

import dr.evolution.io.Importer;
import dr.evolution.io.NexusImporter;
import dr.evolution.io.TreeImporter;
import dr.evolution.tree.Tree;
import dr.evomodelxml.LoggerParser;
import dr.inference.trace.LogFileTraces;
import dr.inference.trace.TraceException;
import dr.stats.DiscreteStatistics;
import dr.util.HeapSort;
import dr.util.TabularData;
import dr.xml.*;

import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.Arrays;

/**
 * @author Joseph Heled
 */
public class EBSPAnalysis extends TabularData {

    private final double[] xPoints;
    private final double[] means;
    private final double[] medians;
    private final double[][] hpdLower;
    private final double[][] hpdHigh;
    private final double[] HPDLevels;
    // each bin covers xPoints[-1]/coalBins.length
    private int[] coalBins;

    private final boolean quantiles;

    EBSPAnalysis(File log, File[] treeFiles, VariableDemographicModel.Type modelType,
                 String firstColumnName, String firstIndicatorColumnName,
                 String rootHeightColumnName, int coalPointBins, double burnIn,
                 double[] inHPDLevels, boolean quantiles, boolean logSpace, boolean mid,
                 int restrictToNchanges)
            throws IOException, Importer.ImportException, TraceException {

        LogFileTraces ltraces = new LogFileTraces(log.getCanonicalPath(), log);
        ltraces.loadTraces();
        ltraces.setBurnIn(0);
        final int runLengthIncludingBurnin = ltraces.getStateCount();

        int intBurnIn = (int) Math.floor(burnIn < 1 ? runLengthIncludingBurnin * burnIn : burnIn);
        final int nStates = runLengthIncludingBurnin - intBurnIn;
        //intBurnIn *= ltraces.getStepSize();
        ltraces.setBurnIn(intBurnIn * ltraces.getStepSize());

        assert ltraces.getStateCount() == nStates;

        this.quantiles = quantiles;
        HPDLevels = (inHPDLevels != null) ? inHPDLevels : new double[]{0.95};

        int populationFirstColumn = -1;
        int indicatorsFirstColumn = -1;
        int rootHeightColumn = -1;

        for (int n = 0; n < ltraces.getTraceCount(); ++n) {
            final String traceName = ltraces.getTraceName(n);
            if (traceName.equals(firstColumnName)) {
                populationFirstColumn = n;
            } else if (traceName.equals(firstIndicatorColumnName)) {
                indicatorsFirstColumn = n;
            } else if (rootHeightColumnName != null && traceName.equals(rootHeightColumnName)) {
                rootHeightColumn = n;
            }
        }

        if (populationFirstColumn < 0 || indicatorsFirstColumn < 0) {
            throw new TraceException("incorrect trace column names: unable to find populations/indicators");
        }

        double binSize = 0;
        if (coalPointBins > 0) {
            if (rootHeightColumn < 0) {
                throw new TraceException("incorrect tree height column");
            }
            double hSum = -0;
            double[] h = new double[1];
            for (int ns = 0; ns < nStates; ++ns) {
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

        for (int k = 0; k < treeFiles.length; ++k) {

            // System.err.println("burnin " + treeFiles[k] + "(" + k + ")");

            treeImporters[k] = new NexusImporter(new FileReader(treeFiles[k]));
            assert intBurnIn > 0;
            for (int z = 0; z < intBurnIn - 1; ++z) {
                treeImporters[k].importNextTree();
            }
            nIndicators += treeImporters[k].importNextTree().getExternalNodeCount() - 1;
        }

        if (isStepWise) {
            nIndicators -= 1;
        }

        final int nXaxisPoints = nIndicators + (isStepWise ? 1 : 0) + 1;
        xPoints = new double[nXaxisPoints];
        Arrays.fill(xPoints, 0.0);

        int nDataPoints = 0;
        VDdemographicFunction[] allDemog = new VDdemographicFunction[nStates];
        {
            double[] indicators = new double[nIndicators];
            double[] pop = new double[nIndicators + 1];
            Tree[] tt = new Tree[treeFiles.length];



            boolean match = true;
            for (int ns = 0; ns < nStates; ++ns) {

                ltraces.getStateValues(ns, indicators, indicatorsFirstColumn);
                ltraces.getStateValues(ns, pop, populationFirstColumn);

                if(match){
                    for (int nt = 0; nt < tt.length; ++nt) {
                        tt[nt] = treeImporters[nt].importNextTree();

                    }
                }
                //Get tree state number
                String name1 = tt[0].getId();
                int state1 = Integer.parseInt(name1.substring(name1.indexOf('_')+1, name1.length()));
                int state2 = state1;

                if (tt.length > 1) {
                    String name2 = tt[1].getId();
                    state2 = Integer.parseInt(name1.substring(name2.indexOf('_')+1, name2.length()));
                }

                if (state1 != state2){     //... can this happen at all?
                    throw new  TraceException("NEXUS tree files have different rates or corrupted!!!!"); //Not too sure what kind of message is appropriate here.


                }else if((ns+intBurnIn)*ltraces.getStepSize() == state1){                   //Check if log state matches tree state
                    match = true;
                    final VDdemographicFunction demoFunction =
                            new VDdemographicFunction(tt, modelType, indicators, pop, logSpace, mid);

                    if (demoFunction.numberOfChanges() == restrictToNchanges) {
                        continue;
                    }

                    double[] xs = demoFunction.allTimePoints();
                    for (int k = 0; k < xs.length; ++k) {
                        xPoints[k + 1] += xs[k];
                    }
                    if (coalPointBins > 0) {
                        for (double x : xs) {
                            coalBins[Math.min((int) (x / binSize), coalBins.length - 1)]++;
                        }
                    }
                    allDemog[nDataPoints] = demoFunction;
                    ++nDataPoints;

                    demoFunction.freeze();
                }else{
                    match = false;
                }

            }

            for (int k = 0; k < xPoints.length; ++k) {
                xPoints[k] /= nStates;
            }
            if(nStates != nDataPoints){                                                     //Warning if log file ant tree files
                // have different rates
                System.err.println("Different Rates is \"main\" and \"tree\" log files");

            }
            if(nDataPoints < 10){                                                           //Warning if number of states is not sufficient
                // enough to do the analysis
                System.err.println("Warning!!! Not Sufficient number of data points");

            }

        }

        double[] popValues = new double[nDataPoints];
        means = new double[nXaxisPoints];
        medians = new double[nXaxisPoints];
        hpdLower = new double[HPDLevels.length][];
        hpdHigh = new double[HPDLevels.length][];

        for (int i = 0; i < HPDLevels.length; ++i) {
            hpdLower[i] = new double[nXaxisPoints];
            hpdHigh[i] = new double[nXaxisPoints];
        }

        for (int nx = 0; nx < xPoints.length; ++nx) {
            final double x = xPoints[nx];

            for (int ns = 0; ns < nDataPoints; ++ns) {
                popValues[ns] = allDemog[ns].getDemographic(x);
            }
            int[] indices = new int[popValues.length];
            HeapSort.sort(popValues, indices);

            means[nx] = DiscreteStatistics.mean(popValues);
            for (int i = 0; i < HPDLevels.length; ++i) {
                if (quantiles) {
                    hpdLower[i][nx] = DiscreteStatistics.quantile((1 - HPDLevels[i]) / 2, popValues, indices);
                    hpdHigh[i][nx] = DiscreteStatistics.quantile((1 + HPDLevels[i]) / 2, popValues, indices);
                } else {
                    final double[] hpd = DiscreteStatistics.HPDInterval(HPDLevels[i], popValues, indices);
                    hpdLower[i][nx] = hpd[0];
                    hpdHigh[i][nx] = hpd[1];
                }
            }
            medians[nx] = DiscreteStatistics.median(popValues, indices);
        }
    }

    private final String[] columnNames = {"time", "mean", "median"};

    public int nColumns() {
        return columnNames.length + 2 * HPDLevels.length + (coalBins != null ? 1 : 0);
    }

    public String columnName(int nColumn) {
        final int fixed = columnNames.length;
        if (nColumn < fixed) {
            return columnNames[nColumn];
        }
        nColumn -= fixed;
        if (nColumn < 2 * HPDLevels.length) {
            final double p = HPDLevels[nColumn / 2];
            final String s = (nColumn % 2 == 0) ? "lower" : "upper";
            return (quantiles ? "cpd " : "hpd ") + s + " " + Math.round(p * 100);
        }
        assert (nColumn - 2 * HPDLevels.length) == 0;
        return "bins";
    }

    public int nRows() {
        return Math.max(xPoints.length, (coalBins != null ? coalBins.length : 0));
    }

    public Object data(int nRow, int nColumn) {
        switch (nColumn) {
            case 0: {
                if (nRow < xPoints.length) {
                    return xPoints[nRow];
                }
                break;
            }
            case 1: {
                if (nRow < means.length) {
                    return means[nRow];
                }
                break;
            }
            case 2: {
                if (nRow < medians.length) {
                    return medians[nRow];
                }
                break;
            }
            default: {
                final int j = nColumn - columnNames.length;
                if (j < 2 * HPDLevels.length) {
                    if (nRow < xPoints.length) {
                        final int k = j / 2;
                        if (0 <= k && k < HPDLevels.length) {
                            if (j % 2 == 0) {
                                return hpdLower[k][nRow];
                            } else {
                                return hpdHigh[k][nRow];
                            }
                        }
                    }
                } else {
                    if (nRow < coalBins.length) {
                        return coalBins[nRow];
                    }
                }
                break;
            }
        }
        return "";
    }

    // should be local to PARSER, but this makes them non-accessible from the outside in Java.

    public static final String VD_ANALYSIS = "VDAnalysis";
    public static final String FILE_NAME = "fileName";
    public static final String BURN_IN = "burnIn";
    public static final String HPD_LEVELS = "Confidencelevels";
    public static final String QUANTILES = "useQuantiles";
    public static final String LOG_SPACE = VariableDemographicModel.LOG_SPACE;
    public static final String USE_MIDDLE = VariableDemographicModel.USE_MIDPOINTS;
    public static final String N_CHANGES = "nChanges";

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
            return xo.getChild(childName).getStringChild(0);
        }

        public String getParserName() {
            return VD_ANALYSIS;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            try {

                // 10% is brun-in default
                final double burnin = xo.getAttribute(BURN_IN, 0.1);
                if (burnin < 0)
                    throw new XMLParseException("burnIn should be either between 0 and 1 or a positive number");

                final double[] hpdLevels = xo.hasAttribute(HPD_LEVELS) ? xo.getDoubleArrayAttribute(HPD_LEVELS) : null;

                final File log = LoggerParser.getFile(getElementText(xo, LOG_FILE_NAME));
                final XMLObject treeFileNames = xo.getChild(TREE_FILE_NAMES);
                final int nTrees = treeFileNames.getChildCount();
                File[] treeFiles = new File[nTrees];
                for (int k = 0; k < nTrees; ++k) {
                    treeFiles[k] = LoggerParser.getFile(((XMLObject) treeFileNames.getChild(k)).getStringChild(0));
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

                final boolean quantiles = xo.getAttribute(QUANTILES, false);
                final boolean logSpace = xo.getAttribute(LOG_SPACE, false);
                final boolean useMid = xo.getAttribute(USE_MIDDLE, false);
                final int onlyNchanges = xo.getAttribute(N_CHANGES, -1);

                return new EBSPAnalysis(log, treeFiles, modelType,
                        populationFirstColumn, indicatorsFirstColumn,
                        rootHeightColumn, nBins,
                        burnin, hpdLevels, quantiles, logSpace, useMid, onlyNchanges);

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
        };
    };
}
