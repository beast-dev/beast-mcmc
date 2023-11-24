/*
 * GaussianMarkovRandomField.java
 *
 * Copyright (c) 2002-2023 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.math.distributions.gp;

import dr.inference.loggers.LogColumn;
import dr.inference.loggers.Loggable;
import dr.inference.loggers.NumberColumn;
import dr.inference.model.*;
import dr.math.MathUtils;
import org.ejml.data.DenseMatrix64F;
import org.ejml.factory.LinearSolverFactory;
import org.ejml.interfaces.linsol.LinearSolver;

import static dr.math.distributions.gp.AdditiveGaussianProcessDistribution.BasisDimension;
import static dr.math.distributions.gp.AdditiveGaussianProcessDistribution.computeAdditiveGramian;

import java.util.ArrayList;
import java.util.List;

/**
 * @author Marc Suchard
 * @author Filippo Monti
 */
public class GaussianProcessPrediction implements Loggable, VariableListener, ModelListener {

    private final AdditiveGaussianProcessDistribution gp;
    private final Parameter realizedValues;
    private final List<DesignMatrix> predictiveDesigns;
    private final int realizedDim;
    private final int predictiveDim;
    private final int order;
    private final double[] prediction;

    private final LinearSolver<DenseMatrix64F> solver;

    private final List<BasisDimension> realizedBases;
    private final List<BasisDimension> predictiveBases;
    private final List<BasisDimension> crossBases;

    private boolean predictionKnown;
    private LogColumn[] columns;

    public GaussianProcessPrediction(AdditiveGaussianProcessDistribution gp,
                                     Parameter realizedValues,
                                     List<DesignMatrix> predictiveDesigns) {
        this.gp = gp;
        this.realizedValues = realizedValues;
        this.predictiveDesigns = predictiveDesigns;
        this.realizedDim = gp.getDimension();
        this.predictiveDim = predictiveDesigns.get(0).getColumnDimension();
        this.order = gp.getOrder();
        this.prediction = new double[realizedDim];

        this.solver = LinearSolverFactory.symmPosDef(realizedDim);

        this.realizedBases = gp.getBases();
        this.predictiveBases = makePredictiveBases(realizedBases, predictiveDesigns);
        this.crossBases = makeCrossBases(realizedBases, predictiveDesigns);

        gp.addModelListener(this);
        realizedValues.addVariableListener(this);

        for (DesignMatrix design : predictiveDesigns) {
            design.addVariableListener(this);
        }

        predictionKnown = false;
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

    private void computePredictions() {
        // Compute: predicition ~ p(f(predictivePoints) | realizedValues)

        DenseMatrix64F realizedGramian = new DenseMatrix64F(realizedDim, realizedDim);
        computeAdditiveGramian(realizedGramian, realizedBases, order); // TODO can get directly from gp

        DenseMatrix64F predictiveGramian = new DenseMatrix64F(predictiveDim, predictiveDim);
        computeAdditiveGramian(predictiveGramian, predictiveBases, order);

        DenseMatrix64F crossGramian = new DenseMatrix64F(realizedDim, predictiveDim);
        computeAdditiveGramian(crossGramian, crossBases, order);

        for (int i = 0; i < predictiveDim; ++i) {
            prediction[i] = MathUtils.nextGaussian();
        }

        // mean: crossGramian %*% inverse(realizedGramian) %*% realizedField
        // variance: predictiveGramian - crossGramian %*% inverse(realizedGramian) %*% t(crossGramian)

        // TODO
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
    public void modelChangedEvent(Model model, Object object, int index) {
        if (model == gp) {
            predictionKnown = false;
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
        } else if (variable instanceof DesignMatrix &&
                predictiveDesigns.contains((DesignMatrix) variable)) {
            throw new IllegalArgumentException("Not yet implemented");
        } else {
            throw new IllegalArgumentException("Unknown variable");
        }
    }
}
