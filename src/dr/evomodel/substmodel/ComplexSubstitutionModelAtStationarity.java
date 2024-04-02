/*
 * ComplexSubstitutionModel.java
 *
 * Copyright (c) 2002-2022 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.evomodel.substmodel;

import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.impl.DenseDoubleMatrix2D;
import cern.colt.matrix.linalg.LUDecomposition;
import dr.evolution.datatype.DataType;
import dr.inference.loggers.LogColumn;
import dr.inference.loggers.NumberColumn;
import dr.inference.model.*;
import dr.util.Citable;
import dr.util.Citation;
import dr.util.CommonCitations;

import java.util.*;

/**
 * @author Jiansi Gao
 * @author Marc Suchard
 */
public class ComplexSubstitutionModelAtStationarity extends ComplexSubstitutionModel {

    public ComplexSubstitutionModelAtStationarity(String name, DataType dataType, Parameter parameter) {
        super(name, dataType, null, parameter);

        stationaryDistribution = new double[stateCount];
        storedStationaryDistribution = new double[stateCount];
        stationaryDistributionKnown = false;
        storedStationaryDistributionKnown = false;
    }

    @Override
    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        stationaryDistributionKnown = false;
        super.handleVariableChangedEvent(variable, index, type);
    }

    @Override
    protected double[] getPi() {
        return getStationaryDistribution();
    }
    
    double[] getStationaryDistribution() {
        if (!stationaryDistributionKnown) {
            computeStationaryDistribution(stationaryDistribution);
            stationaryDistributionKnown = true;
        }
        return stationaryDistribution;
    }

    @Override
    protected void setupQMatrix(double[] rates, double[] pi, double[][] matrix) {
        int i, j, k = 0;
        for (i = 0; i < stateCount; i++) {
            for (j = i + 1; j < stateCount; j++) {
                double thisRate = rates[k++];
                if (thisRate < 0.0) thisRate = 0.0;
                matrix[i][j] = thisRate;
            }
        }
        // Copy lower triangle in column-order form (transposed)
        for (j = 0; j < stateCount; j++) {
            for (i = j + 1; i < stateCount; i++) {
                double thisRate = rates[k++];
                if (thisRate < 0.0) thisRate = 0.0;
                matrix[i][j] = thisRate;
            }
        }
    }

    @Override
    protected void storeState() {
        System.arraycopy(stationaryDistribution, 0, storedStationaryDistribution, 0, stateCount);
        storedStationaryDistributionKnown = stationaryDistributionKnown;

        super.storeState();
    }

    @Override
    protected void restoreState() {
        double[] tmp = stationaryDistribution;
        stationaryDistribution = storedStationaryDistribution;
        storedStationaryDistribution = tmp;
        stationaryDistributionKnown = storedStationaryDistributionKnown;

        super.restoreState();
    }

    private double[] stationaryDistribution;
    private double[] storedStationaryDistribution;
    private boolean stationaryDistributionKnown;
    private boolean storedStationaryDistributionKnown;
}
