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

package dr.inference.operators.factorAnalysis;

import dr.inference.distribution.DistributionLikelihood;
import dr.inference.operators.GibbsOperator;
import dr.inference.operators.SimpleMCMCOperator;
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
 * @author Max R. Tolkoff
 * @author Marc A. Suchard
 */
public class NewLoadingsGibbsOperator extends SimpleMCMCOperator implements GibbsOperator {

    private NormalDistribution workingPrior;
    private final ArrayList<double[][]> precisionArray;
    private final ArrayList<double[]> meanMidArray;
    private final ArrayList<double[]> meanArray;

    private final boolean randomScan;
    private double pathParameter = 1.0;

    private final double priorPrecision;
    private final double priorMean;
    private final double priorPrecisionWorking;

    private final FactorAnalysisOperatorAdaptor adaptor;

    private final ConstrainedSampler constrainedSampler;

    public NewLoadingsGibbsOperator(FactorAnalysisOperatorAdaptor adaptor, DistributionLikelihood prior,
                                    double weight, boolean randomScan, DistributionLikelihood workingPrior,
                                    boolean multiThreaded, int numThreads,
                                    ConstrainedSampler constrainedSampler) {

        setWeight(weight);

        this.adaptor = adaptor;

        NormalDistribution prior1 = (NormalDistribution) prior.getDistribution();
        if (workingPrior != null) {
            this.workingPrior = (NormalDistribution) workingPrior.getDistribution();
        }

        precisionArray = new ArrayList<double[][]>();
        double[][] temp;
        this.randomScan = randomScan;
        this.constrainedSampler = constrainedSampler;

        meanArray = new ArrayList<double[]>();
        meanMidArray = new ArrayList<double[]>();

        // TODO Clean up memory allocation
        double[] tempMean;
        if (!randomScan) {
            for (int i = 0; i < adaptor.getNumberOfFactors(); i++) {
                temp = new double[i + 1][i + 1];
                precisionArray.add(temp);
            }
            for (int i = 0; i < adaptor.getNumberOfFactors(); i++) {
                tempMean = new double[i + 1];
                meanArray.add(tempMean);
            }

            for (int i = 0; i < adaptor.getNumberOfFactors(); i++) {
                tempMean = new double[i + 1];
                meanMidArray.add(tempMean);
            }
        } else {
            for (int i = 0; i < adaptor.getNumberOfFactors(); i++) {
                temp = new double[adaptor.getNumberOfFactors() - i][adaptor.getNumberOfFactors() - i];
                precisionArray.add(temp);
            }
            for (int i = 0; i < adaptor.getNumberOfFactors(); i++) {
                tempMean = new double[adaptor.getNumberOfFactors() - i];
                meanArray.add(tempMean);
            }

            for (int i = 0; i < adaptor.getNumberOfFactors(); i++) {
                tempMean = new double[adaptor.getNumberOfFactors() - i];
                meanMidArray.add(tempMean);
            }
        }

        priorPrecision = 1 / (prior1.getSD() * prior1.getSD());
        priorMean = prior1.getMean();

        if (workingPrior == null) {
            priorPrecisionWorking = priorPrecision;
        } else {
            priorPrecisionWorking = 1 / (this.workingPrior.getSD() * this.workingPrior.getSD());
        }

        if (multiThreaded) {
            for (int i = 0; i < adaptor.getNumberOfTraits(); i++) {
                if (i < adaptor.getNumberOfFactors())
                    drawCallers.add(new DrawCaller(i, new double[i + 1][i + 1], new double[i + 1], new double[i + 1]));
                else
                    drawCallers.add(new DrawCaller(i, new double[adaptor.getNumberOfFactors()][adaptor.getNumberOfFactors()],
                            new double[adaptor.getNumberOfFactors()], new double[adaptor.getNumberOfFactors()]));
            }
            pool = Executors.newFixedThreadPool(numThreads);
        } else {
            pool = null;
        }
    }


    private void getPrecisionOfTruncated(FactorAnalysisOperatorAdaptor adaptor, //MatrixParameterInterface full,
                                         int newRowDimension, int row, double[][] answer) {

        int p = adaptor.getNumberOfTaxa(); //.getColumnDimension();

        for (int i = 0; i < newRowDimension; i++) {
            for (int j = i; j < newRowDimension; j++) {
                double sum = 0;
                for (int k = 0; k < p; k++)
                    sum += adaptor.getFactorValue(i, k) * adaptor.getFactorValue(j, k);
                answer[i][j] = sum * this.adaptor.getColumnPrecision(row); //adaptor.getColumnPrecision().getParameterValue(row, row);
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

        final int p = adaptor.getNumberOfTaxa();

        for (int i = 0; i < newRowDimension; i++) {
            double sum = 0;

            for (int k = 0; k < p; k++) {
                if (adaptor.isNotMissing(dataColumn, k)) {
                    sum += adaptor.getFactorValue(i, k) /*Left.getParameterValue(i, k)*/
                            * adaptor.getDataValue(dataColumn, k); //data.getParameterValue(dataColumn, k);
                }
            }

            sum = sum * adaptor.getColumnPrecision(dataColumn); //adaptor.getColumnPrecision().getParameterValue(dataColumn, dataColumn);
            sum += priorMean * priorPrecision;
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
        int size = adaptor.getNumberOfFactors();
        if (i < size) {
            getPrecisionOfTruncated(adaptor, i + 1, i, answer);
        } else {
            getPrecisionOfTruncated(adaptor, size, i, answer);
        }
    }

    private void getMean(int i, double[][] variance, double[] midMean, double[] mean) {

        int size = adaptor.getNumberOfFactors();
        if (i < size) {
            getTruncatedMean(i + 1, i, variance, midMean, mean);
        } else {
            getTruncatedMean(size, i, variance, midMean, mean);
        }
        for (int j = 0; j < mean.length; j++) {
            mean[j] *= pathParameter;  // TODO Is this missing the working prior component?
        }
    }

    private void drawI(int i, double[][] precision, double[] midMean, double[] mean) {  // TODO Flatten precision

        getPrecision(i, precision);

        double[][] variance = (new SymmetricMatrix(precision)).inverse().toComponents();
        double[][] cholesky = null;
        try {
            cholesky = new CholeskyDecomposition(variance).getL();
        } catch (IllegalDimension illegalDimension) {
            illegalDimension.printStackTrace();
        }

        getMean(i, variance, midMean, mean);

        // TODO Use back-solve to avoid inverting precision first
//                double[] draws = MultivariateNormalDistribution.nextMultivariateNormalViaBackSolvePrecision(
//                         mean, precision);

        double[] draw = MultivariateNormalDistribution.nextMultivariateNormalCholesky(mean, cholesky);

        adaptor.setLoadingsForTraitQuietly(i, draw);

        if (DEBUG) {
            System.err.println("draw: " + new Vector(draws));
        }
    }
    
    @Override
    public String getPerformanceSuggestion() {
        return null;
    }

    @Override
    public String getOperatorName() { return "newLoadingsGibbsOperator"; }

    private static boolean DEBUG = false;

    @Override
    public double doOperation() {

        if (DEBUG) {
            System.err.println("Start doOp");
        }

        // Draw new factors if necessary
        adaptor.drawFactors();

        int size = adaptor.getNumberOfTraits();
        if (adaptor.getNumberOfFactors() != precisionArray.listIterator().next().length) {

            if (DEBUG) {
                System.err.println("!= length");
            }

            precisionArray.clear();
            meanArray.clear();
            meanMidArray.clear();
            double[] tempMean;
            double[][] temp;
            if (!randomScan) {
                for (int i = 0; i < adaptor.getNumberOfFactors(); i++) {
                    temp = new double[i + 1][i + 1];
                    precisionArray.add(temp);
                }
                for (int i = 0; i < adaptor.getNumberOfFactors(); i++) {
                    tempMean = new double[i + 1];
                    meanArray.add(tempMean);
                }

                for (int i = 0; i < adaptor.getNumberOfFactors(); i++) {
                    tempMean = new double[i + 1];
                    meanMidArray.add(tempMean);
                }
            } else {
                for (int i = 0; i < adaptor.getNumberOfFactors(); i++) {
                    temp = new double[adaptor.getNumberOfFactors() - i][adaptor.getNumberOfFactors() - i];
                    precisionArray.add(temp);
                }
                for (int i = 0; i < adaptor.getNumberOfFactors(); i++) {
                    tempMean = new double[adaptor.getNumberOfFactors() - i];
                    meanArray.add(tempMean);
                }

                for (int i = 0; i < adaptor.getNumberOfFactors(); i++) {
                    tempMean = new double[adaptor.getNumberOfFactors() - i];
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
                adaptor.fireLoadingsChanged();
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
                    if (i < adaptor.getNumberOfFactors()) {
                        precision = currentPrecision.next();
                        midMean = currentMidMean.next();
                        mean = currentMean.next();
                    }

                    drawI(i, precision, midMean, mean);
                }
                constrainedSampler.applyConstraint(adaptor);
                adaptor.fireLoadingsChanged();
            } else {
                int i = MathUtils.nextInt(adaptor.getNumberOfTraits());
                ListIterator<double[][]> currentPrecision;
                ListIterator<double[]> currentMidMean;
                ListIterator<double[]> currentMean;
                if (i < adaptor.getNumberOfFactors()) {
                    currentPrecision = precisionArray.listIterator(adaptor.getNumberOfFactors() - i - 1);
                    currentMidMean = meanMidArray.listIterator(adaptor.getNumberOfFactors() - i - 1);
                    currentMean = meanArray.listIterator(adaptor.getNumberOfFactors() - i - 1);
                } else {
                    currentPrecision = precisionArray.listIterator();
                    currentMidMean = meanMidArray.listIterator();
                    currentMean = meanArray.listIterator();
                }
                drawI(i, currentPrecision.next(), currentMidMean.next(), currentMean.next());
                constrainedSampler.applyConstraint(adaptor);
                adaptor.fireLoadingsChanged();
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

        private final static boolean DEBUG_PARALLEL_EVALUATION = false;

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

    public enum ConstrainedSampler {

        NONE("none") {
            @Override
            void applyConstraint(FactorAnalysisOperatorAdaptor adaptor) {
                // Do nothing
            }
        },
        REFLECTION("reflection") {
            @Override
            void applyConstraint(FactorAnalysisOperatorAdaptor adaptor) {
                for (int factor = 0; factor < adaptor.getNumberOfFactors(); ++factor) {
                    adaptor.reflectLoadingsForFactor(factor);
                }
            }
        };

        ConstrainedSampler(String name) {
            this.name = name;
        }

        private String name;

        public String getName() {
            return name;
        }

        public static ConstrainedSampler parse(String name) {
            name = name.toLowerCase();
            for (ConstrainedSampler sampler : ConstrainedSampler.values()) {
                if (name.compareTo(sampler.getName()) == 0) {
                    return sampler;
                }
            }
            throw new IllegalArgumentException("Unknown sampler type");
        }

        abstract void applyConstraint(FactorAnalysisOperatorAdaptor adaptor);
    }
}
