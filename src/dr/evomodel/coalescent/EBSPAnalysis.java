/*
 * EBSPAnalysis.java
 *
 * Copyright (c) 2002-2015 Alexei Drummond, Andrew Rambaut and Marc Suchard
 *
 * This file is part of BEAST.
 * See the NOTICE file distributed with this work for additional
 * information regarding copyright ownership and licensing.
 *
 * BEAST is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 *  BEAST is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with BEAST; if not, write to the
 * Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
 * Boston, MA  02110-1301  USA
 */

package dr.evomodel.coalescent;

import dr.evolution.io.Importer;
import dr.evolution.io.NexusImporter;
import dr.evolution.io.TreeImporter;
import dr.evolution.tree.Tree;
import dr.inference.trace.LogFileTraces;
import dr.inference.trace.TraceException;
import dr.stats.DiscreteStatistics;
import dr.util.HeapSort;
import dr.util.TabularData;

import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
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

    public EBSPAnalysis(File log, File[] treeFiles, VariableDemographicModel.Type modelType,
                        String firstColumnName, String firstIndicatorColumnName,
                        String rootHeightColumnName, int coalPointBins, double burnIn,
                        double[] inHPDLevels, boolean quantiles, boolean logSpace, boolean mid,
                        int restrictToNchanges, PrintWriter allDemoWriter)
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

                if (match) {
                    for (int nt = 0; nt < tt.length; ++nt) {
                        tt[nt] = treeImporters[nt].importNextTree();
                        if( tt[nt] == null ) {
                           throw new TraceException("All NEXUS tree files should contain the same number of states");
                        }
                    }
                }
                //Get tree state number
                final String name1 = tt[0].getId();
                final int state1 = Integer.parseInt(name1.substring(name1.indexOf('_') + 1, name1.length()));

                for (int j = 1; j < tt.length; ++j) {
                    final String name2 = tt[j].getId();
                    int state2 = Integer.parseInt(name1.substring(name2.indexOf('_') + 1, name2.length()));
                    if (state1 != state2) {
                        throw new TraceException("NEXUS tree files have different rates or corrupted!!!!");
                    }
                }

                if ((ns + intBurnIn) * ltraces.getStepSize() == state1) {                   //Check if log state matches tree state
                    match = true;
                    final VDdemographicFunction demoFunction =
                            new VDdemographicFunction(tt, modelType, indicators, pop, logSpace, mid);

                    if (restrictToNchanges >= 0 && demoFunction.numberOfChanges() != restrictToNchanges) {
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
                } else {
                    match = false;
                }
            }

            for (int k = 0; k < xPoints.length; ++k) {
                xPoints[k] /= nStates;
            }

            if (nStates != nDataPoints) {                                                     //Warning if log file and tree files
                // have different rates
                System.err.println("Different Rates is \"main\" and \"tree\" log files");

            }
            if (nDataPoints < 10) {                                                           //Warning if number of states is not sufficient
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

        if( allDemoWriter != null ) {
            for(double xPoint : xPoints) {
                allDemoWriter.print(xPoint);
                allDemoWriter.append(' ');
            }

            for (int ns = 0; ns < nDataPoints; ++ns) {
                allDemoWriter.println();
                for(double xPoint : xPoints) {
                    allDemoWriter.print(allDemog[ns].getDemographic(xPoint));
                    allDemoWriter.append(' ');
                }
            }
            allDemoWriter.close();
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

}
