/*
 * LoadingsGibbsOperator.java
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

package dr.inference.operators;

import dr.inference.distribution.DistributionLikelihood;
import dr.inference.model.LatentFactorModel;
import dr.inference.model.MatrixParameterInterface;
import dr.inference.model.Parameter;
import dr.math.MathUtils;
import dr.math.distributions.MultivariateNormalDistribution;
import dr.math.distributions.NormalDistribution;
import dr.math.matrixAlgebra.*;

import java.util.ArrayList;
import java.util.List;
import java.util.ListIterator;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

/**
 * Author Max R. Tolkoff
 */
public class LoadingsGibbsOperator extends SimpleMCMCOperator implements GibbsOperator {

    private NormalDistribution workingPrior;
    private LatentFactorModel LFM;
    private ArrayList<double[][]> precisionArray;
    private ArrayList<double[]> meanMidArray;
    private ArrayList<double[]> meanArray;
    private boolean randomScan;
    private double pathParameter = 1.0;
    private final Parameter missingIndicator;


    private double priorPrecision;
    private double priorMeanPrecision;
    private double priorPrecisionWorking;

    public LoadingsGibbsOperator(LatentFactorModel LFM, DistributionLikelihood prior, double weight, boolean randomScan, DistributionLikelihood workingPrior, boolean multiThreaded, int numThreads) {
        setWeight(weight);

        NormalDistribution prior1 = (NormalDistribution) prior.getDistribution();
        if (workingPrior != null) {
            this.workingPrior = (NormalDistribution) workingPrior.getDistribution();
        }
        this.LFM = LFM;
        precisionArray = new ArrayList<double[][]>();
        double[][] temp;
        this.randomScan = randomScan;


        meanArray = new ArrayList<double[]>();
        meanMidArray = new ArrayList<double[]>();
        double[] tempMean;
        if (!randomScan) {
            for (int i = 0; i < LFM.getFactorDimension(); i++) {
                temp = new double[i + 1][i + 1];
                precisionArray.add(temp);
            }
            for (int i = 0; i < LFM.getFactorDimension(); i++) {
                tempMean = new double[i + 1];
                meanArray.add(tempMean);
            }

            for (int i = 0; i < LFM.getFactorDimension(); i++) {
                tempMean = new double[i + 1];
                meanMidArray.add(tempMean);
            }
        } else {
            for (int i = 0; i < LFM.getFactorDimension(); i++) {
                temp = new double[LFM.getFactorDimension() - i][LFM.getFactorDimension() - i];
                precisionArray.add(temp);
            }
            for (int i = 0; i < LFM.getFactorDimension(); i++) {
                tempMean = new double[LFM.getFactorDimension() - i];
                meanArray.add(tempMean);
            }

            for (int i = 0; i < LFM.getFactorDimension(); i++) {
                tempMean = new double[LFM.getFactorDimension() - i];
                meanMidArray.add(tempMean);
            }
        }

//            vectorProductAnswer=new MatrixParameter[LFM.getLoadings().getRowDimension()];
//            for (int i = 0; i <vectorProductAnswer.length ; i++) {
//                vectorProductAnswer[i]=new MatrixParameter(null);
//                vectorProductAnswer[i].setDimensions(i+1, 1);
//            }

//        priorMeanVector=new MatrixParameter[LFM.getLoadings().getRowDimension()];
//            for (int i = 0; i <priorMeanVector.length ; i++) {
//                priorMeanVector[i]=new MatrixParameter(null, i+1, 1, this.prior.getMean()/(this.prior.getSD()*this.prior.getSD()));
//
//
//            }
        priorPrecision = 1 / (prior1.getSD() * prior1.getSD());
        priorMeanPrecision = prior1.getMean() * priorPrecision;

//        double priorMeanPrecisionWorking;
        if (workingPrior == null) {
//            priorMeanPrecisionWorking = priorMeanPrecision;
            priorPrecisionWorking = priorPrecision;
        } else {
            priorPrecisionWorking = 1 / (this.workingPrior.getSD() * this.workingPrior.getSD());
//            priorMeanPrecisionWorking = this.workingPrior.getMean() * priorPrecisionWorking;
        }

        if (multiThreaded) {
            for (int i = 0; i < LFM.getLoadings().getRowDimension(); i++) {
                if (i < LFM.getFactorDimension())
                    drawCallers.add(new DrawCaller(i, new double[i + 1][i + 1], new double[i + 1], new double[i + 1]));
                else
                    drawCallers.add(new DrawCaller(i, new double[LFM.getFactorDimension()][LFM.getFactorDimension()], new double[LFM.getFactorDimension()], new double[LFM.getFactorDimension()]));
            }

//                    Integer.parseInt(System.getProperty("thread.count"));
            pool = Executors.newFixedThreadPool(numThreads);
        } else {
            pool = null;
        }
        missingIndicator = LFM.getMissingIndicator();
    }


    private void getPrecisionOfTruncated(MatrixParameterInterface full, int newRowDimension, int row, double[][] answer) {

//        MatrixParameter answer=new MatrixParameter(null);
//        answer.setDimensions(this.getRowDimension(), Right.getRowDimension());
//        System.out.println(answer.getRowDimension());
//        System.out.println(answer.getColumnDimension());

        int p = full.getColumnDimension();
        for (int i = 0; i < newRowDimension; i++) {
            for (int j = i; j < newRowDimension; j++) {
                double sum = 0;
                for (int k = 0; k < p; k++)
                    sum += full.getParameterValue(i, k) * full.getParameterValue(j, k);
                answer[i][j] = sum * LFM.getColumnPrecision().getParameterValue(row, row);
                if (i == j) {
                    answer[i][j] = answer[i][j] * pathParameter + getAdjustedPriorPrecision();
                } else {
                    answer[i][j] *= pathParameter;
                    answer[j][i] = answer[i][j];
                }
            }
        }
    }


    private void getTruncatedMean(int newRowDimension, int dataColumn, double[][] variance, double[] midMean, double[] mean) {

//        MatrixParameter answer=new MatrixParameter(null);
//        answer.setDimensions(this.getRowDimension(), Right.getRowDimension());
//        System.out.println(answer.getRowDimension());
//        System.out.println(answer.getColumnDimension());

        MatrixParameterInterface data = LFM.getScaledData();
        MatrixParameterInterface Left = LFM.getFactors();

        final int p = data.getColumnDimension();

        for (int i = 0; i < newRowDimension; i++) {
            double sum = 0;

            for (int k = 0; k < p; k++) {
                if (missingIndicator == null || missingIndicator.getParameterValue(k * LFM.getScaledData().getRowDimension() + dataColumn) != 1)
                    sum += Left.getParameterValue(i, k) * data.getParameterValue(dataColumn, k);
            }

            sum = sum * LFM.getColumnPrecision().getParameterValue(dataColumn, dataColumn);
            sum += priorMeanPrecision;
            midMean[i] = sum;
        }

        for (int i = 0; i < newRowDimension; i++) {
            double sum = 0;
            for (int k = 0; k < newRowDimension; k++)
                sum += variance[i][k] * midMean[k];
            mean[i] = sum;
        }
    }

    private void getPrecision(int i, double[][] answer) {
        int size = LFM.getFactorDimension();
        if (i < size) {
            getPrecisionOfTruncated(LFM.getFactors(), i + 1, i, answer);
        } else {
            getPrecisionOfTruncated(LFM.getFactors(), size, i, answer);
        }
    }

    private void getMean(int i, double[][] variance, double[] midMean, double[] mean) {
//        Matrix factors=null;
        int size = LFM.getFactorDimension();
//        double[] scaledDataColumn=LFM.getScaledData().getRowValues(i);
//        Vector dataColumn=null;
//        Vector priorVector=null;
//        Vector temp=null;
//        Matrix data=new Matrix(LFM.getScaledData().getParameterAsMatrix());
        if (i < size) {
            getTruncatedMean(i + 1, i, variance, midMean, mean);
//            dataColumn=new Vector(data.toComponents()[i]);
//            try {
//                answer=precision.inverse().product(new Matrix(priorMeanVector[i].add(vectorProductAnswer[i]).getParameterAsMatrix()));
//            } catch (IllegalDimension illegalDimension) {
//                illegalDimension.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
//            }
        } else {
            getTruncatedMean(size, i, variance, midMean, mean);
//            dataColumn=new Vector(data.toComponents()[i]);
//            try {
//                answer=precision.inverse().product(new Matrix(priorMeanVector[size-1].add(vectorProductAnswer[size-1]).getParameterAsMatrix()));
//            } catch (IllegalDimension illegalDimension) {
//                illegalDimension.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
//            }
        }
        for (int j = 0; j < mean.length; j++) { //TODO implement for generic prior
            mean[j] *= pathParameter;
        }

    }

    private void copy(int i, double[] random) {
        MatrixParameterInterface changing = LFM.getLoadings();
        for (int j = 0; j < random.length; j++) {
            changing.setParameterValueQuietly(i, j, random[j]);
        }
    }

    private void drawI(int i, double[][] precision, double[] midMean, double[] mean) {
        double[] draws;
        double[][] variance;
        double[][] cholesky = null;
        getPrecision(i, precision);
        variance = (new SymmetricMatrix(precision)).inverse().toComponents();

        try {
            cholesky = new CholeskyDecomposition(variance).getL();
        } catch (IllegalDimension illegalDimension) {
            illegalDimension.printStackTrace();
        }

        getMean(i, variance, midMean, mean);

        draws = MultivariateNormalDistribution.nextMultivariateNormalCholesky(mean, cholesky);
//    if(i<draws.length)
//
//    {
//        while (draws[i] < 0) {
//            draws = MultivariateNormalDistribution.nextMultivariateNormalCholesky(mean, cholesky);
//        }
//    }
        if (i < draws.length) {
            //if (draws[i] > 0) { // TODO implement as option
            copy(i, draws);
            //}
        } else {
            copy(i, draws);
        }

        if (DEBUG) {
            System.err.println("draw: " + new Vector(draws));
        }

//       copy(i, draws);

    }
    
    @Override
    public String getPerformanceSuggestion() {
        return null;
    }

    @Override
    public String getOperatorName() {
        return "loadingsGibbsOperator";
    }

    private static boolean DEBUG = false;

    @Override
    public double doOperation() {

        if (DEBUG) {
            System.err.println("Start doOp");
        }

        int size = LFM.getLoadings().getRowDimension();
        if (LFM.getFactorDimension() != precisionArray.listIterator().next().length) {

            if (DEBUG) {
                System.err.println("!= length");
            }

            precisionArray.clear();
            meanArray.clear();
            meanMidArray.clear();
            double[] tempMean;
            double[][] temp;
            if (!randomScan) {
                for (int i = 0; i < LFM.getFactorDimension(); i++) {
                    temp = new double[i + 1][i + 1];
                    precisionArray.add(temp);
                }
                for (int i = 0; i < LFM.getFactorDimension(); i++) {
                    tempMean = new double[i + 1];
                    meanArray.add(tempMean);
                }

                for (int i = 0; i < LFM.getFactorDimension(); i++) {
                    tempMean = new double[i + 1];
                    meanMidArray.add(tempMean);
                }
            } else {
                for (int i = 0; i < LFM.getFactorDimension(); i++) {
                    temp = new double[LFM.getFactorDimension() - i][LFM.getFactorDimension() - i];
                    precisionArray.add(temp);
                }
                for (int i = 0; i < LFM.getFactorDimension(); i++) {
                    tempMean = new double[LFM.getFactorDimension() - i];
                    meanArray.add(tempMean);
                }

                for (int i = 0; i < LFM.getFactorDimension(); i++) {
                    tempMean = new double[LFM.getFactorDimension() - i];
                    meanMidArray.add(tempMean);
                }
            }
        }

        if (pool != null) {

            if (DEBUG) {
                System.err.println("!= poll");
            }

            try {
                pool.invokeAll(drawCallers);
                LFM.getLoadings().fireParameterChangedEvent();
            } catch (InterruptedException e) {
                e.printStackTrace();
            }
        } else {

            if (DEBUG) {
                System.err.println("inner");
            }

            if (!randomScan) {
                ListIterator<double[][]> currentPrecision = precisionArray.listIterator();
                ListIterator<double[]> currentMidMean = meanMidArray.listIterator();
                ListIterator<double[]> currentMean = meanArray.listIterator();
                double[][] precision = null;
                double[] midMean = null;
                double[] mean = null;
                for (int i = 0; i < size; i++) {
                    if (i < LFM.getFactorDimension()) {
                        precision = currentPrecision.next();
                        midMean = currentMidMean.next();
                        mean = currentMean.next();
                    }

                    drawI(i, precision, midMean, mean);
                }
                LFM.getLoadings().fireParameterChangedEvent();
            } else {
                int i = MathUtils.nextInt(LFM.getLoadings().getRowDimension());
                ListIterator<double[][]> currentPrecision;
                ListIterator<double[]> currentMidMean;
                ListIterator<double[]> currentMean;
                if (i < LFM.getFactorDimension()) {
                    currentPrecision = precisionArray.listIterator(LFM.getFactorDimension() - i - 1);
                    currentMidMean = meanMidArray.listIterator(LFM.getFactorDimension() - i - 1);
                    currentMean = meanArray.listIterator(LFM.getFactorDimension() - i - 1);
                } else {
                    currentPrecision = precisionArray.listIterator();
                    currentMidMean = meanMidArray.listIterator();
                    currentMean = meanArray.listIterator();
                }
                drawI(i, currentPrecision.next(), currentMidMean.next(), currentMean.next());
                LFM.getLoadings().fireParameterChangedEvent(i, null);
//            LFM.getLoadings().fireParameterChangedEvent();
            }

        }

        if (DEBUG) {
            for (double[] m : meanArray) {
                System.err.println(new Vector(m));
            }

            for (double[] m : meanMidArray) {
                System.err.println(new Vector(m));
            }

            for (double[][] p : precisionArray) {
                System.err.println(new Matrix(p));
            }
//            System.err.println(new Matrix(precisionArray));
        }

        if (DEBUG) {
            System.err.println("End doOp");
        }

        return 0;
    }

    public void setPathParameter(double beta) {
        pathParameter = beta;
    }

    private double getAdjustedPriorPrecision() {
        return priorPrecision * pathParameter + (1 - pathParameter) * priorPrecisionWorking;
    }


    class DrawCaller implements Callable<Double> {

        int i;
        double[][] precision;
        double[] midMean;
        double[] mean;

        DrawCaller(int i, double[][] precision, double[] midMean, double[] mean) {
            this.i = i;
            this.precision = precision;
            this.midMean = midMean;
            this.mean = mean;
        }

        private final boolean DEBUG_PARALLEL_EVALUATION = false;

        public Double call() throws Exception {

            if (DEBUG_PARALLEL_EVALUATION) {
                System.err.print("Invoking thread #" + i + " for " + ": ");
            }
            drawI(i, precision, midMean, mean);
            return null;
        }

    }

    private final List<Callable<Double>> drawCallers = new ArrayList<Callable<Double>>();

    private final ExecutorService pool;
}
