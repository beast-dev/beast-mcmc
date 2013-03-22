package dr.evomodel.coalescent;

import dr.evolution.io.Importer;
import dr.inference.model.Parameter;

//import dr.evolution.io.NexusImporter;
//import dr.evolution.io.TreeImporter;
//import dr.evolution.tree.Tree;
import dr.inference.operators.CoercionMode;
import dr.inference.trace.LogFileTraces;
import dr.inference.trace.TraceException;
import dr.evomodel.coalescent.operators.GaussianProcessSkytrackBlockUpdateOperator;
import dr.stats.DiscreteStatistics;
import dr.util.FileHelpers;
import dr.util.HeapSort;
import dr.util.TabularData;
import no.uib.cipr.matrix.*;


import java.io.*;
//import java.io.PrintWriter;
import java.util.Arrays;

/**
 * @author Joseph Heled
 */
public class GPSkytrackAnalysis extends TabularData {

    private final double[] xPoints;
    private final double[] means;
    private final double[] medians;
    private final double[] hpdLower;
    private final double[] hpdHigh;
    private final double [][] gValues;
    private final double [][] tValues;
    private final double [][] newGvalues;
    private final double [][] popValues;


//    final File gvalues = FileHelpers.getFile("gvalues.txt");
//    final File locations = FileHelpers.getFile("locations.txt");

    //    private final double[] HPDLevels;
    private Parameter numGridpoints;
    // each bin covers xPoints[-1]/coalBins.length
//    private int[] coalBins;

//    private final boolean quantiles;

//    GaussianProcessSkytrackLikelihood gpLikelihood = (GaussianProcessSkytrackLikelihood) xo.getChild(GaussianProcessSkytrackLikelihood.class);
//    return new GaussianProcessSkytrackBlockUpdateOperator(gpLikelihood, weight, mode, scaleFactor,
//                                                          maxIterations, stopValue);

    public GPSkytrackAnalysis(File log,  double burnIn, Parameter numGridPoints) throws IOException, Importer.ImportException, TraceException {
        GaussianProcessSkytrackBlockUpdateOperator GPOperator=new GaussianProcessSkytrackBlockUpdateOperator();
        this.numGridpoints=numGridPoints;
        LogFileTraces ltraces = new LogFileTraces(log.getCanonicalPath(), log);
        ltraces.loadTraces();

        ltraces.setBurnIn(0);
        final int runLengthIncludingBurnin = ltraces.getStateCount();

        int intBurnIn = (int) Math.floor(burnIn < 1 ? runLengthIncludingBurnin * burnIn : burnIn);
        final int nStates = runLengthIncludingBurnin - intBurnIn;
//        System.err.println("runl"+runLengthIncludingBurnin+"burnin");
        //intBurnIn *= ltraces.getStepSize();
        ltraces.setBurnIn(intBurnIn * ltraces.getStepSize());
        assert ltraces.getStateCount() == nStates;

//        this.quantiles = quantiles;
//        HPDLevels = new double[]{0.95};

        xPoints = new double[(int) numGridPoints.getParameterValue(0)+1];
        means = new double[(int) numGridPoints.getParameterValue(0)+1];
        medians = new double[(int) numGridPoints.getParameterValue(0)+1];
        hpdHigh = new double[(int) numGridPoints.getParameterValue(0)+1];
        hpdLower = new double[(int) numGridPoints.getParameterValue(0)+1];

        int numbPointsColumn = -1;
        int lambdaColumn = -1;
        int precColumn = -1;
        int tmrcaColumn=-1;



        for (int n = 0; n < ltraces.getTraceCount(); ++n) {
            final String traceName = ltraces.getTraceName(n);
            if (traceName.equals("skyride.points")) {
               numbPointsColumn = n;
            } else if (traceName.equals("skyride.lambda_bound")) {
                lambdaColumn = n;
            } else if (traceName.equals("skyride.precision")) {
                precColumn = n;
            } else if (traceName.equals("skyride.tmrca")) {
                tmrcaColumn = n;
            }

        }

        if (numbPointsColumn < 0 || lambdaColumn < 0 || precColumn<0 || tmrcaColumn<0) {
            throw new TraceException("incorrect trace column names: unable to find lambda_bound, points, tmrca or precision");
        }

//        TODO: Check if it is ok to define the grid from 0 to max(TMRCA) always
        double binSize = 0;
//            double hSum = -0;


            int [] numPoints = new int[nStates];
            double[] lambda = new double[nStates];
            double[] kappa= new double[nStates];
            double tmrca=0;
//            double binSize=0;
            double tempTmrca=0.0;
            int maxpts=0;
            for (int ns = 0; ns < nStates; ++ns) {
                lambda[ns]= (Double) ltraces.getTrace(lambdaColumn).getValue(ns);
                numPoints[ns]=(int)Math.round((Double) ltraces.getTrace(numbPointsColumn).getValue(ns));
                kappa[ns]=(Double) ltraces.getTrace(precColumn).getValue(ns);
                tempTmrca=(Double) ltraces.getTrace(tmrcaColumn).getValue(ns);
                if (tempTmrca>tmrca){tmrca=tempTmrca;}
                if (numPoints[ns]>maxpts) {maxpts=numPoints[ns];}
            }

            binSize = tmrca / numGridPoints.getParameterValue(0);
            xPoints[0]=0.0;
            for (int np=1;np<xPoints.length;np++){
             xPoints[np]=xPoints[np-1]+binSize;
            }

            gValues=new double[nStates][];
            tValues=new double[nStates][];
            newGvalues=new double[nStates][];
            popValues=new double[nStates][];
            readChain(gValues,"gvalues.txt");
            readChain(tValues,"locations.txt");

           for (int j=0;j<nStates;j++){
//               newGvalues[j]=new double[numPoints[j]];
               newGvalues[j]=GPOperator.getGPvaluesS(tValues[j], gValues[j], xPoints, kappa[j]);
//               popValues[j]=(1+Math.exp(-newGvalues[j].))/lambda[j];
               }


//            System.err.println(gvalTraces.getTrace(0).getValues(0,5));


//             for (int j=0;j<nStates)
//        GaussianProcessSkytrackBlockUpdateOperator.getGPvaluesS(currentChangePoints, currentPopSize, xPoints,currentPrecision);


//
//        int nDataPoints = 0;
//        VDdemographicFunction[] allDemog = new VDdemographicFunction[nStates];
//        {
//            double[] indicators = new double[nIndicators];
//            double[] pop = new double[nIndicators + 1];
//            Tree[] tt = new Tree[treeFiles.length];
//
//            boolean match = true;
//            for (int ns = 0; ns < nStates; ++ns) {
//
//                ltraces.getStateValues(ns, indicators, indicatorsFirstColumn);
//                ltraces.getStateValues(ns, pop, populationFirstColumn);
//
//                if (match) {
//                    for (int nt = 0; nt < tt.length; ++nt) {
//                        tt[nt] = treeImporters[nt].importNextTree();
//                        if( tt[nt] == null ) {
//                           throw new TraceException("All NEXUS tree files should contain the same number of states");
//                        }
//                    }
//                }
//                //Get tree state number
//                final String name1 = tt[0].getId();
//                final int state1 = Integer.parseInt(name1.substring(name1.indexOf('_') + 1, name1.length()));
//
//                for (int j = 1; j < tt.length; ++j) {
//                    final String name2 = tt[j].getId();
//                    int state2 = Integer.parseInt(name1.substring(name2.indexOf('_') + 1, name2.length()));
//                    if (state1 != state2) {
//                        throw new TraceException("NEXUS tree files have different rates or corrupted!!!!");
//                    }
//                }
//
//                if ((ns + intBurnIn) * ltraces.getStepSize() == state1) {                   //Check if log state matches tree state
//                    match = true;
//                    final VDdemographicFunction demoFunction =
//                            new VDdemographicFunction(tt, modelType, indicators, pop, logSpace, mid);
//
//                    if (restrictToNchanges >= 0 && demoFunction.numberOfChanges() != restrictToNchanges) {
//                        continue;
//                    }
//
//                    double[] xs = demoFunction.allTimePoints();
//                    for (int k = 0; k < xs.length; ++k) {
//                        xPoints[k + 1] += xs[k];
//                    }
//                    if (coalPointBins > 0) {
//                        for (double x : xs) {
//                            coalBins[Math.min((int) (x / binSize), coalBins.length - 1)]++;
//                        }
//                    }
//                    allDemog[nDataPoints] = demoFunction;
//                    ++nDataPoints;
//
//                    demoFunction.freeze();
//                } else {
//                    match = false;
//                }
//            }
//
//            for (int k = 0; k < xPoints.length; ++k) {
//                xPoints[k] /= nStates;
//            }
//
//            if (nStates != nDataPoints) {                                                     //Warning if log file and tree files
//                // have different rates
//                System.err.println("Different Rates is \"main\" and \"tree\" log files");
//
//            }
//            if (nDataPoints < 10) {                                                           //Warning if number of states is not sufficient
//                // enough to do the analysis
//                System.err.println("Warning!!! Not Sufficient number of data points");
//            }
//        }
//
//        double[] popValues = new double[nDataPoints];
//        means = new double[nXaxisPoints];
//        medians = new double[nXaxisPoints];
//        hpdLower = new double[HPDLevels.length][];
//        hpdHigh = new double[HPDLevels.length][];
//
//        for (int i = 0; i < HPDLevels.length; ++i) {
//            hpdLower[i] = new double[nXaxisPoints];
//            hpdHigh[i] = new double[nXaxisPoints];
//        }
//
//        for (int nx = 0; nx < xPoints.length; ++nx) {
//            final double x = xPoints[nx];
//
//            for (int ns = 0; ns < nDataPoints; ++ns) {
//                popValues[ns] = allDemog[ns].getDemographic(x);
//            }
//            int[] indices = new int[popValues.length];
//            HeapSort.sort(popValues, indices);
//
//            means[nx] = DiscreteStatistics.mean(popValues);
//            for (int i = 0; i < HPDLevels.length; ++i) {
//                if (quantiles) {
//                    hpdLower[i][nx] = DiscreteStatistics.quantile((1 - HPDLevels[i]) / 2, popValues, indices);
//                    hpdHigh[i][nx] = DiscreteStatistics.quantile((1 + HPDLevels[i]) / 2, popValues, indices);
//                } else {
//                    final double[] hpd = DiscreteStatistics.HPDInterval(HPDLevels[i], popValues, indices);
//                    hpdLower[i][nx] = hpd[0];
//                    hpdHigh[i][nx] = hpd[1];
//                }
//            }
//            medians[nx] = DiscreteStatistics.median(popValues, indices);
//        }
//
//        if( allDemoWriter != null ) {
//            for(double xPoint : xPoints) {
//                allDemoWriter.print(xPoint);
//                allDemoWriter.append(' ');
//            }
//
//            for (int ns = 0; ns < nDataPoints; ++ns) {
//                allDemoWriter.println();
//                for(double xPoint : xPoints) {
//                    allDemoWriter.print(allDemog[ns].getDemographic(xPoint));
//                    allDemoWriter.append(' ');
//                }
//            }
//            allDemoWriter.close();
//        }
    }

    public void readChain(double [][] current,String fileName){
        try {
            BufferedReader br = new BufferedReader(new FileReader(fileName));

            String line=null;
            int i=0;
//            System.err.println("will read line1");
            while ((line=br.readLine())!=null){
                String[] parts=line.split(" ");
//                System.err.println(i+"with cols:"+parts.length);
                current[i]=new double[parts.length];
                    for (int j=0;j<parts.length;j++){
                        current[i][j]=Double.parseDouble(parts[j]);
                    }
                i++;
        }   br.close();
        }
        catch(java.io.IOException ioe){
            System.err.println("IOException:"+ ioe.getMessage());
        }
    }



    private final String[] columnNames = {"time", "mean", "median","upper","lower"};

    public int nColumns() {
        return 5;
    }

    public String columnName(int nColumn) {
//        final int fixed = columnNames.length;
//        if (nColumn < fixed) {
            return columnNames[nColumn];
        }
//        nColumn -= fixed;
//        if (nColumn < 2 * HPDLevels.length) {
//            final double p = HPDLevels[nColumn / 2];
//            final String s = (nColumn % 2 == 0) ? "lower" : "upper";
//            return (quantiles ? "cpd " : "hpd ") + s + " " + Math.round(p * 100);
//        }
//        assert (nColumn - 2 * HPDLevels.length) == 0;
//        return "bins";
//    }

    public int nRows() {
        return (int) numGridpoints.getParameterValue(0)+1;
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
            case 3: {
                if (nRow < hpdLower.length) {
                    return hpdLower[nRow];
                }
                break;
            }
            case 4: {
                if (nRow < hpdHigh.length) {
                    return hpdHigh[nRow];
                }
                break;
            }
//            default: {
//                final int j = nColumn - columnNames.length;
//                if (j < 2 * HPDLevels.length) {
//                    if (nRow < xPoints.length) {
//                        final int k = j / 2;
//                        if (0 <= k && k < HPDLevels.length) {
//                            if (j % 2 == 0) {
//                                return hpdLower[k][nRow];
//                            } else {
//                                return hpdHigh[k][nRow];
//                            }
//                        }
//                    }
//                } else {
//                    if (nRow < coalBins.length) {
//                        return coalBins[nRow];
//                    }
//                }
//                break;
//            }
        }
        return "";
    }

}
