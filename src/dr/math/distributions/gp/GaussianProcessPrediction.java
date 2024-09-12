/*
 * GaussianProcessPrediction.java
 *
 * Copyright © 2002-2024 the BEAST Development Team
 * http://beast.community/about
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
 *
 */

package dr.math.distributions.gp;

import dr.inference.loggers.LogColumn;
import dr.inference.loggers.Loggable;
import dr.inference.loggers.NumberColumn;
import dr.inference.model.*;
import dr.math.MathUtils;
import dr.xml.Reportable;
import org.ejml.alg.dense.decomposition.chol.CholeskyDecompositionCommon_D64;
import org.ejml.data.DenseMatrix64F;
import org.ejml.factory.LinearSolverFactory;
import org.ejml.interfaces.linsol.LinearSolver;
import org.ejml.ops.CommonOps;

import static dr.math.distributions.gp.AdditiveGaussianProcessDistribution.BasisDimension;
import static dr.math.distributions.gp.AdditiveGaussianProcessDistribution.computeAdditiveGramian;

import java.util.ArrayList;
import java.util.List;

/**
 * @author Marc Suchard
 * @author Filippo Monti
 */
public class GaussianProcessPrediction implements Reportable, Loggable, VariableListener, ModelListener {

    private final AdditiveGaussianProcessDistribution gp;
    private final Parameter realizedValues;
    private final List<DesignMatrix> predictiveDesigns;
    private final int realizedDim;
    private final int predictiveDim;
    private final Parameter orderVariance;
    private final double[] prediction;
    private final double[] mean;
    private final DenseMatrix64F variance;

    private final DenseMatrix64F crossGramian;
    private DenseMatrix64F realizedPrecision;

    private final LinearSolver<DenseMatrix64F> solver;

    private final List<BasisDimension> predictiveBases;
    private final List<BasisDimension> crossBases;

    private final DenseMatrix64F crossRealized;

    private boolean predictionKnown;
    private boolean meanKnown;
    private boolean varianceKnown;
    private boolean crossRealizedKnown; // crossRealized = crossGramian%*%inverse(realizedGramian) TODO I am actually using the realizedPrecision directly
    private boolean crossGramianKnown;
    private LogColumn[] columns;

    public GaussianProcessPrediction(AdditiveGaussianProcessDistribution gp,
                                     Parameter realizedValues,
                                     List<DesignMatrix> predictiveDesigns) {
        this.gp = gp;
        this.realizedValues = realizedValues;
        this.predictiveDesigns = predictiveDesigns;
        this.realizedDim = gp.getDimension();
        this.predictiveDim = predictiveDesigns.get(0).getRowDimension(); // CHANGED TO ROW ??? (BEFORE IT WAS ColumnDimension)
        this.orderVariance = gp.getOrderVariance();
        this.crossGramian = new DenseMatrix64F(predictiveDim, realizedDim);
        this.realizedPrecision = new DenseMatrix64F(realizedDim, realizedDim);
        this.crossRealized = new DenseMatrix64F(predictiveDim, realizedDim);

        this.mean = new double[predictiveDim];
        this.variance = new DenseMatrix64F(predictiveDim,predictiveDim);
        this.prediction = new double[predictiveDim]; // CHANGED FROM  realizedDim

        this.solver = LinearSolverFactory.symmPosDef(realizedDim);

        List<BasisDimension> realizedBases = gp.getBases();
        this.predictiveBases = makePredictiveBases(realizedBases, predictiveDesigns);
        this.crossBases = makeCrossBases(realizedBases, predictiveDesigns);

        gp.addModelListener(this);
        realizedValues.addVariableListener(this);

        for (DesignMatrix design : predictiveDesigns) {
            design.addVariableListener(this);
        }

        predictionKnown = false;
        meanKnown = false;
        varianceKnown = false;
        crossRealizedKnown = false;
        crossGramianKnown = false;
    }

    private static List<BasisDimension> makeCrossBases(List<BasisDimension> originalBases,
                                                            List<DesignMatrix> predictiveDesigns) {
        List<BasisDimension> result = new ArrayList<>();

        for (int i = 0; i < originalBases.size(); ++i) {
            BasisDimension originalBasis = originalBases.get(i);
            BasisDimension newBasis = new BasisDimension(
                    originalBasis.getKernel(), predictiveDesigns.get(i), originalBasis.getDesignMatrix1());

            result.add(newBasis);
        }

        return result;
    }

    private static List<BasisDimension> makePredictiveBases(List<BasisDimension> originalBases,
                                                            List<DesignMatrix> predictiveDesigns) {
        List<BasisDimension> result = new ArrayList<>();

        for (int i = 0; i < originalBases.size(); ++i) {
            BasisDimension originalBasis = originalBases.get(i);
            BasisDimension newBasis = new BasisDimension(
                    originalBasis.getKernel(), predictiveDesigns.get(i), predictiveDesigns.get(i));

            result.add(newBasis);
        }

        return result;
    }

    // CHECK nextMultivariateNormalViaBackSolvePrecision

    private void computePredictions() {
        // Compute: prediction ~ p(f(predictivePoints) | realizedValues)
        computeMean();
        computeVariance();

        // TODO here I am overwriting the variance with the lowerTriangular matrix. I do not want that
        if (!solver.setA(variance)) {
            throw new RuntimeException("Unable to decompose matrix");
        }
        DenseMatrix64F lowerTriangularVar =  ((CholeskyDecompositionCommon_D64) solver.getDecomposition()).getT();

        double[] standardGaussian = new double[predictiveDim];
        // multivariate standard normal generator
        for (int i = 0; i < predictiveDim; ++i) {
            standardGaussian[i] = MathUtils.nextGaussian();
        }
        // computing the prediction
        for(int i = 0; i < predictiveDim; i++) {
            prediction[i] = 0; // re-setting the prediction
            for(int j = 0; j < predictiveDim; j++) {
                prediction[i] += lowerTriangularVar.get(i,j) * standardGaussian[j];
            }
            prediction[i] += mean[i];
        }
    }

    private double getPrediction(int index) {
        if (!predictionKnown) {
            computePredictions();
            predictionKnown = true;
        }

        if (index == predictiveDim - 1) {
            predictionKnown = false;
        }
        return prediction[index];
    }

    private void computeCrossGramian() {
        if (!crossGramianKnown) {
            computeAdditiveGramian(crossGramian, crossBases, orderVariance);
            crossGramianKnown = true;
        }
    }

    private void computeCrossRealized() {
        if (!crossRealizedKnown) {
            computeCrossGramian();

            realizedPrecision = gp.getPrecisionAsMatrix();
            CommonOps.mult(crossGramian, realizedPrecision, crossRealized);
            crossRealizedKnown = true;
        }
    }

    private void computeMean() {
        if (!meanKnown) {
            computeCrossRealized();
            for(int i = 0; i < predictiveDim; i += 1) {
                mean[i] = 0;
                for(int j = 0; j < realizedDim; j += 1) {
                    mean[i] += crossRealized.get(i,j) * realizedValues.getParameterValue(j);
                }
            }
            meanKnown = true;
        }
    }

    private void computeVariance() {
        if (!varianceKnown) {
            computeCrossRealized();

            DenseMatrix64F predictiveGramian = new DenseMatrix64F(predictiveDim, predictiveDim);
            computeAdditiveGramian(predictiveGramian, predictiveBases, orderVariance);

            DenseMatrix64F predictiveGramianCorrection = new DenseMatrix64F(predictiveDim, predictiveDim);
            for (int i = 0; i < crossRealized.numRows; i++) {
                for (int j = 0; j < crossGramian.numRows; j++) {
                    double sum = 0.0;
                    for (int k = 0; k < crossRealized.numCols; k++) {
                        sum += crossRealized.get(i, k) * crossGramian.get(j, k); // multiplying for the transpose
                    }
                    predictiveGramianCorrection.set(i, j, sum);
                }
            }
            // variance: predictiveGramian - [crossGramian %*% inverse(realizedGramian) %*% t(crossGramian)]
            CommonOps.subtract(predictiveGramian, predictiveGramianCorrection, variance);
        }
    }

    private double[] getMean() {
        computeMean();
        return mean;
    }

    private double[] getVariance() {
        // TODO decide if I want to compute variance or directly the lowerTriangular Cholensky matrix
        computeVariance();
        return variance.getData();
    }

    private LogColumn[] createColumns() {
        LogColumn[] columns = new LogColumn[predictiveDim];
        for (int i = 0; i < predictiveDim; ++i) {
            final String name = "prediction" + (i + 1);
            final int index = i;
            columns[i] = new NumberColumn(name) {
                @Override
                public double getDoubleValue() {
                    return getPrediction(index);
                }
            };
        }

        return columns;
    }

    @Override
    public LogColumn[] getColumns() {
        if (columns == null) {
            columns = createColumns();
        }
        return columns;
    }

    @Override
    public String getReport() {
        StringBuilder sb = new StringBuilder();
        sb.append("mean:");
        for (double value : getMean()) {
            sb.append(" ").append(value);
        }
        sb.append("\n");
        sb.append("variance:");
        for (double value : getVariance()) {
            sb.append(" ").append(value);
        }
        sb.append("\n");
        sb.append("prediction:");
        for (int i = 0; i < predictiveDim;  i++) {
            sb.append(" ").append(getPrediction(i));
        }
        return sb.toString();
    }

    @Override
    public void modelChangedEvent(Model model, Object object, int index) {
        if (model == gp) {
            predictionKnown = false;
            meanKnown = false;
            varianceKnown = false;
            crossRealizedKnown = false;
            crossGramianKnown = false;
        } else {
            throw new IllegalArgumentException("Unknown model");
        }
    }

    @Override
    public void modelRestored(Model model) { predictionKnown = false; }

    @Override
    public void variableChangedEvent(Variable variable, int index, Variable.ChangeType type) {
        if (variable == realizedValues) {
            predictionKnown = false;
            meanKnown = false;
        } else if (variable instanceof DesignMatrix &&
                predictiveDesigns.contains((DesignMatrix) variable)) {
            throw new IllegalArgumentException("Not yet implemented");
        } else {
            throw new IllegalArgumentException("Unknown variable");
        }
    }
}