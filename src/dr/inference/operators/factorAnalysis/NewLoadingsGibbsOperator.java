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

import dr.evomodel.treedatalikelihood.continuous.HashedMissingArray;
import dr.inference.distribution.DistributionLikelihood;
import dr.inference.distribution.NormalDistributionModel;
import dr.inference.distribution.NormalStatisticsProvider;
import dr.inference.operators.GibbsOperator;
import dr.inference.operators.SimpleMCMCOperator;
import dr.math.MathUtils;
import dr.math.distributions.MultivariateNormalDistribution;
import dr.math.matrixAlgebra.*;
import dr.xml.Reportable;
import org.ejml.data.DenseMatrix64F;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

/**
 * @author Max R. Tolkoff
 * @author Marc A. Suchard
 * @author Gabriel Hassler
 */
public class NewLoadingsGibbsOperator extends SimpleMCMCOperator implements GibbsOperator, Reportable {

    private final boolean useInnerProductCache;
    private Map<HashedMissingArray, DenseMatrix64F> precisionMatrixMap = new HashMap<>();


    private NormalDistributionModel workingPrior;
    private final ArrayList<double[][]> precisionArray;
    private final ArrayList<double[]> meanMidArray;
    private final ArrayList<double[]> meanArray;

    private final boolean randomScan;
    private double pathParameter = 1.0;

    private final NormalStatisticsProvider prior;
    private final double priorPrecisionWorking;

    private final FactorAnalysisOperatorAdaptor adaptor;

    private final ConstrainedSampler constrainedSampler;
    private final ColumnDimProvider columnDimProvider;

    private final double[][] observedIndicators;

    public NewLoadingsGibbsOperator(FactorAnalysisOperatorAdaptor adaptor, NormalStatisticsProvider prior,
                                    double weight, boolean randomScan, DistributionLikelihood workingPrior,
                                    boolean multiThreaded, int numThreads,
                                    ConstrainedSampler constrainedSampler,
                                    ColumnDimProvider columnDimProvider,
                                    CacheProvider cacheProvider
    ) {

        setWeight(weight);

        this.adaptor = adaptor;

        this.prior = prior;

        if (workingPrior != null) {
            this.workingPrior = (NormalDistributionModel) workingPrior.getDistribution();
        }

        precisionArray = new ArrayList<double[][]>();
        meanMidArray = new ArrayList<double[]>();
        meanArray = new ArrayList<double[]>();

        this.randomScan = randomScan;
        this.constrainedSampler = constrainedSampler;
        this.columnDimProvider = columnDimProvider;


        if (workingPrior == null) {
            priorPrecisionWorking = getPrecision(prior, 0); //TODO: Do I need to let priorPrecisionWorking vary with dimension?
        } else {
            priorPrecisionWorking = 1 / (this.workingPrior.getStdev() * this.workingPrior.getStdev());
        }

        if (multiThreaded) {
            for (int i = 0; i < adaptor.getNumberOfTraits(); i++) {
                int dim = columnDimProvider.getColumnDim(i, adaptor.getNumberOfFactors());
                drawCallers.add(new DrawCaller(i, new double[dim][dim], new double[dim], new double[dim]));
            }
            pool = Executors.newFixedThreadPool(numThreads);
        } else {
            pool = null;
            columnDimProvider.allocateStorage(precisionArray, meanMidArray, meanArray, adaptor.getNumberOfFactors());
        }

        this.useInnerProductCache = cacheProvider.useCache();

        if (useInnerProductCache) {
            if (multiThreaded && numThreads > 1) {
                throw new IllegalArgumentException("Cannot currently parallelize cached precisions");
            }

            observedIndicators = setupObservedIndicators();

        } else {
            observedIndicators = null;
        }
    }

    public enum CacheProvider {
        USE_CACHE {
            @Override
            boolean useCache() {
                return true;
            }

        },
        NO_CACHE {
            @Override
            boolean useCache() {
                return false;
            }
        };

        abstract boolean useCache();

    }

    private double getPrecision(NormalStatisticsProvider provider, int dim) {
        double sd = provider.getNormalSD(dim);
        return 1.0 / (sd * sd);
    }

    private double[][] setupObservedIndicators() {
        double[][] obsInds = new double[adaptor.getNumberOfTraits()][adaptor.getNumberOfTaxa()];

        for (int trait = 0; trait < adaptor.getNumberOfTraits(); trait++) {
            for (int taxon = 0; taxon < adaptor.getNumberOfTaxa(); taxon++) {

                if (adaptor.isNotMissing(trait, taxon)) {
                    obsInds[trait][taxon] = 1.0;
                }
            }
        }


        return obsInds;
    }


    private void getPrecisionOfTruncated(FactorAnalysisOperatorAdaptor adaptor, //MatrixParameterInterface full,
                                         int newRowDimension, int row, double[][] answer) {

        DenseMatrix64F hashedPrecision = null;
        HashedMissingArray observedArray = null;

        if (useInnerProductCache) {
            double[] observed = observedIndicators[row];
            observedArray = new HashedMissingArray(observed);
            hashedPrecision = precisionMatrixMap.get(observedArray);
        }

        if (!useInnerProductCache || hashedPrecision == null) {

            int p = adaptor.getNumberOfTaxa(); //.getColumnDimension();

            for (int i = 0; i < newRowDimension; i++) {
                for (int j = i; j < newRowDimension; j++) {
                    double sum = 0;
                    for (int k = 0; k < p; k++) {
                        if (adaptor.isNotMissing(row, k)) {
                            sum += adaptor.getFactorValue(i, k) * adaptor.getFactorValue(j, k);
                        }
                    }

                    answer[i][j] = sum;
                    if (i != j) {
                        answer[j][i] = sum;
                    }
                }
            }

            if (useInnerProductCache) {
                precisionMatrixMap.put(observedArray, new DenseMatrix64F(answer));
            }
        } else {
            for (int i = 0; i < newRowDimension; i++) {
                System.arraycopy(hashedPrecision.getData(), i * newRowDimension,
                        answer[i], 0, newRowDimension);
            }

        }

        for (int i = 0; i < newRowDimension; i++) {
            for (int j = i; j < newRowDimension; j++) {
                answer[i][j] *= this.adaptor.getColumnPrecision(row); //adaptor.getColumnPrecision().getParameterValue(row, row);
                if (i == j) {
                    answer[i][j] = answer[i][j] * pathParameter +
                            getAdjustedPriorPrecision(adaptor.getNumberOfTraits() * i + row);
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

            int priorDim = adaptor.getNumberOfTraits() * i + dataColumn;

            sum = sum * adaptor.getColumnPrecision(dataColumn); //adaptor.getColumnPrecision().getParameterValue(dataColumn, dataColumn);
            sum += prior.getNormalMean(priorDim) * getPrecision(prior, priorDim);
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
        getPrecisionOfTruncated(adaptor, columnDimProvider.getColumnDim(i, size), i, answer);
    }

    private void getMean(int i, double[][] variance, double[] midMean, double[] mean) {

        int size = adaptor.getNumberOfFactors();
        getTruncatedMean(columnDimProvider.getColumnDim(i, size), i, variance, midMean, mean);

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

        double[] draw = MultivariateNormalDistribution.nextMultivariateNormalCholesky(mean, cholesky);

        adaptor.setLoadingsForTraitQuietly(i, draw);

        if (DEBUG) {
            System.err.println("draw: " + new Vector(draw));
        }
    }

    private void drawI(int i) {
        int arrayInd = columnDimProvider.getArrayIndex(i, adaptor.getNumberOfFactors());
        drawI(i, precisionArray.get(arrayInd), meanMidArray.get(arrayInd), meanArray.get(arrayInd));
    }

//    @Override
//    public String getPerformanceSuggestion() {
//        return null;
//    }

    @Override
    public String getOperatorName() {
        return "newLoadingsGibbsOperator";
    }

    private static boolean DEBUG = false;


    @Override
    public double doOperation() {

        if (DEBUG) {
            System.err.println("Start doOp");
        }

        // Draw new factors if necessary
        adaptor.drawFactors();

        int size = adaptor.getNumberOfTraits();

        if (useInnerProductCache) {
            precisionMatrixMap.clear();
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
                for (int i = 0; i < size; i++) {
                    drawI(i);
                }

            } else {
                int i = MathUtils.nextInt(adaptor.getNumberOfTraits());
                drawI(i);

            }

            constrainedSampler.applyConstraint(adaptor);
            adaptor.fireLoadingsChanged();

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

    private double getAdjustedPriorPrecision(int dim) {
        return getPrecision(prior, dim) * pathParameter + (1 - pathParameter) * priorPrecisionWorking;
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

    public enum ColumnDimProvider {


        NONE("none") {
            @Override
            int getColumnDim(int colIndex, int nRows) {
                return nRows;
            }

            @Override
            int getArrayIndex(int colIndex, int nRows) {
                return 0;
            }

            @Override
            void allocateStorage(ArrayList<double[][]> precisionArray, ArrayList<double[]> midMeanArray,
                                 ArrayList<double[]> meanArray, int nRows) {

                precisionArray.add(new double[nRows][nRows]);
                midMeanArray.add(new double[nRows]);
                meanArray.add(new double[nRows]);

            }
        },

        UPPER_TRIANGULAR("upperTriangular") {
            @Override
            int getColumnDim(int colIndex, int nRows) {
                return Math.min(colIndex + 1, nRows);
            }

            @Override
            int getArrayIndex(int colIndex, int nRows) {
                return Math.min(colIndex, nRows - 1);
            }

            @Override
            void allocateStorage(ArrayList<double[][]> precisionArray, ArrayList<double[]> midMeanArray,
                                 ArrayList<double[]> meanArray, int nRows) {

                for (int i = 1; i <= nRows; i++) {
                    precisionArray.add(new double[i][i]);
                    midMeanArray.add(new double[i]);
                    meanArray.add(new double[i]);
                }

            }
        };


        abstract int getColumnDim(int colIndex, int nRows);

        abstract int getArrayIndex(int colIndex, int nRows);

        abstract void allocateStorage(ArrayList<double[][]> precisionArray, ArrayList<double[]> midMeanArray,
                                      ArrayList<double[]> meanArray, int nRows);


        private String name;

        ColumnDimProvider(String name) {
            this.name = name;
        }

        public String getName() {
            return name;
        }

        public static ColumnDimProvider parse(String name) {
            name = name.toLowerCase();
            for (ColumnDimProvider dimProvider : ColumnDimProvider.values()) {
                if (name.compareTo(dimProvider.getName().toLowerCase()) == 0) {
                    return dimProvider;
                }
            }
            throw new IllegalArgumentException("Unknown dimension provider type");
        }

    }

    @Override
    public String getReport() {
        int repeats = 1000000;
        int nFac = adaptor.getNumberOfFactors();
        int nTraits = adaptor.getNumberOfTraits();


        int dimLoadings = nTraits * nFac;
        double[] loadMean = new double[dimLoadings];
        double[][] loadCov = new double[dimLoadings][dimLoadings];

        double[] originalLoadings = new double[dimLoadings];
        for (int i = 0; i < dimLoadings; i++) {
            originalLoadings[i] = adaptor.getLoadingsValue(i);
        }


        for (int rep = 0; rep < repeats; rep++) {
            doOperation();
            for (int i = 0; i < dimLoadings; i++) {
                loadMean[i] += adaptor.getLoadingsValue(i);
                for (int j = i; j < dimLoadings; j++) {
                    loadCov[i][j] += adaptor.getLoadingsValue(i) * adaptor.getLoadingsValue(j);
                }
            }
            adaptor.fireLoadingsChanged();
        }

        restoreLoadings(originalLoadings);
        adaptor.fireLoadingsChanged();

        for (int i = 0; i < dimLoadings; i++) {
            loadMean[i] /= repeats;
            for (int j = i; j < dimLoadings; j++) {
                loadCov[i][j] /= repeats;
            }
        }


        for (int i = 0; i < dimLoadings; i++) {
            for (int j = i; j < dimLoadings; j++) {
                loadCov[i][j] = loadCov[i][j] - loadMean[i] * loadMean[j];
                loadCov[j][i] = loadCov[i][j];
            }
        }

        StringBuilder sb = new StringBuilder();
        sb.append(getOperatorName() + "Report:\n");
        sb.append("Loadings mean:\n");
        sb.append(new Vector(loadMean));
        sb.append("\n\n");
        sb.append("Loadings covariance:\n");
        sb.append(new Matrix(loadCov));
        sb.append("\n\n");

        return sb.toString();
    }

    private void restoreLoadings(double[] originalLoadings) {
        int nTraits = adaptor.getNumberOfTraits();
        int nFac = adaptor.getNumberOfFactors();

        double[] buffer = new double[nFac];

        for (int i = 0; i < nTraits; i++) {

            for (int j = 0; j < nFac; j++) {
                buffer[j] = originalLoadings[j * nTraits + i];
            }

            adaptor.setLoadingsForTraitQuietly(i, buffer);
        }
    }

}

