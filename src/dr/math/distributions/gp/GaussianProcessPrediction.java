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

/**
 * @author Marc Suchard
 * @author Filippo Monti
 */
public class GaussianProcessPrediction implements Loggable, VariableListener, ModelListener {

    private final GaussianProcessDistribution gp;
    private final Parameter realizedValues;
    private final Parameter predictivePoints;
    private final int dim;
    private final double[] prediction;

    private boolean predictionKnown;
    private LogColumn[] columns;

    public GaussianProcessPrediction(GaussianProcessDistribution gp,
                                     Parameter realizedValues,
                                     Parameter predictivePoints) {
        this.gp = gp;
        this.realizedValues = realizedValues;
        this.predictivePoints = predictivePoints;
        this.dim = gp.getDimension();
        this.prediction = new double[dim];

        predictionKnown = false;
    }

    private void computePredictions() {
        // TODO
        // Compute: predicition ~ p(f(predictivePoints) | realizedValues)
    }

    private double getPrediction(int index) {
        if (!predictionKnown) {
            computePredictions();
            predictionKnown = true;
        }
        return prediction[index];
    }

    private LogColumn[] createColumns() {
        LogColumn[] columns = new LogColumn[dim];
        for (int i = 0; i < dim; ++i) {
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
        if (variable == realizedValues || variable == predictivePoints) {
            predictionKnown = false;
        } else {
            throw new IllegalArgumentException("Unknown variable");
        }
    }
}
