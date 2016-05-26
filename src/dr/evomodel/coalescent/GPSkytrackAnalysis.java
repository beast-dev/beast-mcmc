/*
 * GPSkytrackAnalysis.java
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
import dr.inference.model.Parameter;

//import dr.evolution.io.NexusImporter;
//import dr.evolution.io.TreeImporter;
//import dr.evolution.tree.Tree;
//import dr.inference.operators.CoercionMode;
//import dr.inference.trace.AbstractTraceList;
import dr.inference.trace.LogFileTraces;
import dr.inference.trace.TraceException;
import dr.evomodel.coalescent.operators.GaussianProcessSkytrackBlockUpdateOperator;
//import dr.inference.trace.TraceFactory;
import dr.stats.DiscreteStatistics;
//import dr.util.FileHelpers;
//import dr.util.HeapSort;
import dr.util.TabularData;
//import no.uib.cipr.matrix.SymmTridiagEVD;
//import no.uib.cipr.matrix.*;


import java.io.*;
//import java.io.PrintWriter;
import java.util.Arrays;
import java.util.StringTokenizer;

/**
 * @author Joseph Heled
 */
public class GPSkytrackAnalysis extends TabularData {
//    TabularData

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
    private Parameter numGridPoints;
    // each bin covers xPoints[-1]/coalBins.length
//    private int[] coalBins;

//    private final boolean quantiles;

//    GaussianProcessSkytrackLikelihood gpLikelihood = (GaussianProcessSkytrackLikelihood) xo.getChild(GaussianProcessSkytrackLikelihood.class);
//    return new GaussianProcessSkytrackBlockUpdateOperator(gpLikelihood, weight, mode, scaleFactor,
//                                                          maxIterations, stopValue);

//    TODO: Error in loadTraces() because a String {..} is being converted to "real/double"
//    To make my life more miserable I will not use logFileTraces class and do it by hand


    public GPSkytrackAnalysis(File log,  double burnIn, Parameter numGridPoints) throws IOException, Importer.ImportException, TraceException {
        GaussianProcessSkytrackBlockUpdateOperator GPOperator=new GaussianProcessSkytrackBlockUpdateOperator();
        this.numGridPoints=numGridPoints;
        LogFileTraces ltraces = new LogFileTraces(log.getCanonicalPath(), log);

//        ltraces.changeTraceType(1, TraceFactory.TraceType.STRING);

        ltraces.loadTraces();
//        System.exit(-1);

        ltraces.setBurnIn(0);
        final int runLengthIncludingBurnin = ltraces.getStateCount();

        int intBurnIn = (int) Math.floor(burnIn < 1 ? runLengthIncludingBurnin * burnIn : burnIn);
        final int nStates = runLengthIncludingBurnin - intBurnIn;
        ltraces.setBurnIn(intBurnIn * ltraces.getStepSize());
        assert ltraces.getStateCount() == nStates;


        xPoints = new double[(int) numGridPoints.getParameterValue(0)+1];
        means = new double[(int) numGridPoints.getParameterValue(0)+1];
        medians = new double[(int) numGridPoints.getParameterValue(0)+1];
        hpdHigh = new double[(int) numGridPoints.getParameterValue(0)+1];
        hpdLower = new double[(int) numGridPoints.getParameterValue(0)+1];

        int numbPointsColumn = -1;
        int gvaluesColumn=-1;
        int xvaluesColumn=-1;
        int lambdaColumn = -1;
        int precColumn = -1;
        int tmrcaColumn=-1;



        for (int n = 0; n < ltraces.getTraceCount(); ++n) {
            final String traceName = ltraces.getTraceName(n);
            System.err.println(traceName);
            if (traceName.equals("skyride.points")) {
               numbPointsColumn = n;
            } else if (traceName.equals("skyride.lambda_bound")) {
                lambdaColumn = n;
            } else if (traceName.equals("skyride.precision")) {
                precColumn = n;
            } else if (traceName.equals("skyride.tmrca")) {
                tmrcaColumn = n;
            } else if (traceName.equals("changePoints")){
                xvaluesColumn=n;
            } else if (traceName.equals("Gvalues")){
                gvaluesColumn=n;
            }


        }
//        System.err.println("columns"+tmrcaColumn+" tmrca"+xvaluesColumn+" and"+gvaluesColumn);
        if (numbPointsColumn < 0 || lambdaColumn < 0 || precColumn<0 || tmrcaColumn<0 || xvaluesColumn<0 || gvaluesColumn<0) {
            throw new TraceException("incorrect trace column names: unable to find correct columns for summary");
        }

//        TODO: Check if it is ok to define the grid from 0 to max(TMRCA) always
        double binSize = 0;
//            double hSum = -0;

//                           System.err.println("states"+nStates);
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
//                System.err.println(tempTmrca);
                System.exit(-1);
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
            popValues=new double[(int) numGridPoints.getParameterValue(0)+1][];
            readChain(gValues,"gvalues.txt");
            readChain(tValues,"locations.txt");
          for (int i=0;i<=numGridPoints.getParameterValue(0);i++){
            popValues[i]=new double[nStates-1] ;
           }
//
           for (int j=0;j<nStates-1;j++){
//               newGvalues[j]=new double[numPoints[j]];

               newGvalues[j]=GPOperator.getGPvaluesS(tValues[j], gValues[j], xPoints, kappa[j]);
//               popValues[j]=new double[nStates];
                for (int i=0;i<=numGridPoints.getParameterValue(0);i++){
                    popValues[i][j]=(1+Math.exp(-newGvalues[j][i]))/lambda[j];
                }
               }
//
////
//        hpdLower = new double[HPDLevels.length][];
//        hpdHigh = new double[HPDLevels.length][];
//

        for (int nx = 0; nx < xPoints.length; ++nx) {
            means[nx] = DiscreteStatistics.mean(popValues[nx]);
            medians[nx]=DiscreteStatistics.median(popValues[nx]);
            hpdLower[nx]=DiscreteStatistics.quantile(0.025,popValues[nx]);
            hpdHigh[nx]=DiscreteStatistics.quantile(0.975,popValues[nx]);

        }



//
//          for (int i = 0; i < HPDLevels.length; ++i) {
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



    private final String[] columnNames = {"time", "mean", "median","lower","upper"};

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
        return (int) numGridPoints.getParameterValue(0)+1;
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
