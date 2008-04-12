package dr.inference.trace;

import dr.xml.*;
import dr.evomodel.tree.TreeTrace;
import dr.evomodel.coalescent.VariableDemographicModel;
import dr.evomodel.coalescent.VDdemographicFunction;
import dr.evolution.io.Importer;
import dr.evolution.tree.Tree;
import dr.stats.DiscreteStatistics;
import dr.util.HeapSort;
import dr.util.TabularData;
import dr.inference.loggers.MCLogger;

import java.io.*;
import java.util.Arrays;

/**
 @author Joseph Heled
 */
public class VariableDemographicAnalysis extends TabularData {

    private double[] xPoints;
    private double[] means;
    private double[] medians ;
    private double[][] hpdLower;
    private double[][] hpdHigh;
    private double[] HPDLevels;
    private boolean quantiles;

    VariableDemographicAnalysis(File log, File[] treeFiles, VariableDemographicModel.Type modelType,
                                String firstColumnName, String firstIndicatorColumnName, double burnIn,
                                double[] inHPDLevels, boolean quantiles)
            throws IOException, Importer.ImportException, TraceException {

        LogFileTraces ltraces = new LogFileTraces(log.getCanonicalPath(), log);
        ltraces.loadTraces();
        ltraces.setBurnIn(0);
        final int runLengthIncludingBurnin = ltraces.getStateCount();

        int intBurnIn =  (int)Math.floor(burnIn < 1 ? runLengthIncludingBurnin * burnIn : burnIn);
        final int nStates = runLengthIncludingBurnin - intBurnIn;
        intBurnIn *= ltraces.getStepSize();
        ltraces.setBurnIn(intBurnIn);

        assert ltraces.getStateCount() == nStates;

        this.quantiles = quantiles;
        HPDLevels = (inHPDLevels != null ) ? inHPDLevels :  new double[]{0.95};

        int populationFirstColumn = -1;
        int indicatorsFirstColumn = -1;

        for(int n = 0; n < ltraces.getTraceCount(); ++n) {
            final String traceName = ltraces.getTraceName(n);
            if( traceName.equals(firstColumnName) ) {
                populationFirstColumn = n;
            } else if( traceName.equals(firstIndicatorColumnName) ) {
                indicatorsFirstColumn = n;
            }
        }

        if( populationFirstColumn < 0 || indicatorsFirstColumn < 0 ) {
            throw new TraceException("incorrect trace column names: unable to find populations/indicators");
        }

        TreeTrace[] traces = new TreeTrace[treeFiles.length];
        final boolean isStepWise = modelType == VariableDemographicModel.Type.STEPWISE;

        int nIndicators = 0;

        for(int k = 0; k < treeFiles.length; ++k) {

            final TreeTrace treeTrace = TreeTrace.loadTreeTrace(new FileReader(treeFiles[k]));

            nIndicators += treeTrace.getTree(0,0).getExternalNodeCount() - 1;

            final int kthRunLength = treeTrace.getTreeCount(0);

            if( runLengthIncludingBurnin != kthRunLength ) {
                throw new  IOException("non matching runs: " + runLengthIncludingBurnin + " != " +  kthRunLength) ; // FIXME another type
            }
            traces[k] = treeTrace;
        }

        if( isStepWise ) {
            nIndicators -= 1;
        }

//        double heightLimit;
//        {
//            double[] allHeights = new double[traces.length * nStates];
//            int nh = 0;
//            for (final TreeTrace treeTrace : traces) {
//                for (int n = 0; n < nStates; ++n) {
//                    final Tree tree = treeTrace.getTree(n, intBurnIn);
//                    final double height = tree.getNodeHeight(tree.getRoot());
//                    allHeights[nh] = height;
//                    ++nh;
//                }
//            }
//            heightLimit = dr.stats.DiscreteStatistics.quantile(xaxisratio, allHeights);
//        }

        final int nXaxisPoints = nIndicators + (isStepWise ? 1 : 0) + 1;
        xPoints = new double[nXaxisPoints];
        Arrays.fill(xPoints, 0.0);

        VDdemographicFunction[] allDemog = new VDdemographicFunction[nStates];
        {
            double[] indicators = new double[nIndicators];
            double[] pop = new double[nIndicators+1];
            Tree[] tt = new Tree[treeFiles.length];
            //double[] popFactors = null;

            for(int ns = 0; ns < nStates; ++ns) {
                ltraces.getStateValues(ns, indicators,  indicatorsFirstColumn);
                ltraces.getStateValues(ns, pop, populationFirstColumn);

                for(int nt = 0; nt < tt.length; ++nt) {
                    tt[nt] = traces[nt].getTree(ns, intBurnIn);
                }
                final VDdemographicFunction demoFunction = new VDdemographicFunction(tt, modelType /*, popFactors*/, indicators, pop);
                double[] x = demoFunction.allTimePoints();
                for(int k = 0; k < x.length; ++k) {
                    assert x[k] >= 0 : " " + k + " " + x[k];
                    xPoints[k+1] += x[k];
                }
                allDemog[ns] = demoFunction;
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
            medians[nx] =  DiscreteStatistics.median(popValues, indices);
//            for(int k = 0; k < indices.length; ++k) {
//                System.out.print(popValues[indices[k]]);  System.out.print(",");
//            }
        }
    }

    private String[] columnNames = {"time", "mean", "median"};
          //, "hpd lower", "hpd upper"};
    public int nColumns() {
        return columnNames.length + 2*HPDLevels.length;
    }


    public String columnName(int nColumn) {
        final int fixed = columnNames.length;
        if( nColumn < fixed) {
            return columnNames[nColumn];
        }
        nColumn -= fixed;
        double p = HPDLevels[nColumn/2];
        String s = ( nColumn % 2 == 0 ) ? "lower" : "upper";
        return (quantiles ? "cpd " : "hpd ") + s + " " + Math.round(p*100);
    }

    public int nRows() {
        return xPoints.length;
    }

    public Object data(int nRow, int nColumn) {
        switch( nColumn ) {
            case 0:
            {
                return xPoints[nRow];
            }
            case 1:
            {
                return means[nRow];
            }
            case 2:
            {
                return medians[nRow];
            }
            default:
            {
                final int j = nColumn - columnNames.length;
                final int k = j /2;
                if( 0 <= k && k < HPDLevels.length ) {
                    if( j % 2 == 0 ) {
                        return hpdLower[k][nRow];
                    } else {
                        return hpdHigh[k][nRow];
                    }
                }
                break;
            }
        }
        return null;
    }

    public static dr.xml.XMLObjectParser PARSER = new dr.xml.AbstractXMLObjectParser() {
        public static final String VD_ANALYSIS = "VDAnalysis";
        public static final String FILE_NAME = "fileName";
        public static final String BURN_IN = "burnIn";
        public static final String HPD_LEVELS = "Confidencelevels";
        public static final String QUANTILES = "useQuantiles";

        public static final String TREE_LOG = "treeOfLoci";

        public static final String LOG_FILE_NAME = "logFileName";
        public static final String TREE_FILE_NAMES = "treeFileNames";
        public static final String MODEL_TYPE = "populationModelType";
        public static final String POPULATION_FIRST_COLUMN = "populationFirstColumn";
        public static final String INDICATORS_FIRST_COLUMN = "indicatorsFirstColumn";


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

                boolean quantiles = xo.hasAttribute(QUANTILES) && xo.getBooleanAttribute(QUANTILES);
                return new VariableDemographicAnalysis(log, treeFiles, modelType, populationFirstColumn,
                        indicatorsFirstColumn, burnin, hpdLevels, quantiles);

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
            return "reconstruct population graph from variable dimension run.";
        }

        public Class getReturnType() {
            return VariableDemographicAnalysis.class;
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
                AttributeRule.newDoubleRule(BURN_IN, true, "The number of states (not sampled states, but" +
                        " actual states) that are discarded from the beginning of the trace before doing the analysis"),
                AttributeRule.newDoubleArrayRule(HPD_LEVELS, true),
                AttributeRule.newBooleanRule(QUANTILES, true),
                
                new ElementRule(LOG_FILE_NAME, String.class, "The name of a BEAST log file"),
                new ElementRule(TREE_FILE_NAMES,
                        new XMLSyntaxRule[] {
                             new ElementRule(TREE_LOG, String.class, "The name of a BEAST trees log file", 1, Integer.MAX_VALUE)
                        } ),
                new ElementRule(MODEL_TYPE, String.class, "population model type (stepwise, linear ..."),
                new ElementRule(POPULATION_FIRST_COLUMN, String.class, "Name of first column of population size"),
                new ElementRule(INDICATORS_FIRST_COLUMN, String.class, "Name of first column of population indicators"),
        };
    };
}
