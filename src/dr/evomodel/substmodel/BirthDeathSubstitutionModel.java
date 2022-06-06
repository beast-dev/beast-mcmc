/*
 * BirthDeathSubstitutionModel.java
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

import dr.evolution.datatype.DataType;
import dr.inference.model.Bounds;
import dr.inference.model.Parameter;

/**
 * @author Marc Suchard
 */
public class BirthDeathSubstitutionModel extends ComplexSubstitutionModel {

    private final Parameter birthParameter;
    private final Parameter deathParameter;

    private final BirthDeathParameterization parameterization;

    public BirthDeathSubstitutionModel(String name,
                                       Parameter birthParameter,
                                       Parameter deathParameter,
                                       DataType dataType) {
        super(name, dataType, null, null);

        this.birthParameter = birthParameter;
        this.deathParameter = deathParameter;

        addVariable(birthParameter);
        addVariable(deathParameter);

        this.parameterization = BirthDeathParameterization.LINEAR;
        this.freqModel = setupEquilibriumModel();

        checkDataType(dataType);
    }

    @Override
    protected void frequenciesChanged() {
        throw new RuntimeException("Frequencies are fixed to birth-death equilibrium distribution");
    }

    @Override
    protected void setupRelativeRates(double[] rates) {
        // Do nothing
    }

    @Override
    protected void setupQMatrix(double[] rates, double[] pi, double[][] matrix) {
        final double perCapitaBirthRate = birthParameter.getParameterValue(0);
        final double perCapitaDeathRate = deathParameter.getParameterValue(0);

        for (int i = 1; i < stateCount; i++) {
            matrix[i - 1][i] = parameterization.birthRate(i,  perCapitaBirthRate);
            matrix[i][i - 1] = parameterization.deathRate(i + 1,  perCapitaDeathRate);
        }
    }

    @Override
    protected double getNormalizationValue(double[][] matrix, double[] pi) {
        return 1.0;
    }

    @Override
    public double getLogLikelihood() {
        return 0.0;
    }
    
    private FrequencyModel setupEquilibriumModel() {
        Parameter equilibrium = new Parameter.Proxy("equilibrium", stateCount) {

            @Override
            public double getParameterValue(int dim) {
                return 1.0 / stateCount; // TODO Fix!
            }

            @Override
            public void setParameterValue(int dim, double value) {
                throw new IllegalArgumentException("Cannot set equilibrium");
            }

            @Override
            public void setParameterValueQuietly(int dim, double value) {
                throw new IllegalArgumentException("Cannot set equilibrium");
            }

            @Override
            public void setParameterValueNotifyChangedAll(int dim, double value) {
                throw new IllegalArgumentException("Cannot set equilibrium");
            }

            @Override
            public void addBounds(Bounds<Double> bounds) {
                // Do nothing
            }
        };

        return new FrequencyModel(dataType, equilibrium);
    }

    private static void checkDataType(DataType dataType) {
        for (int i = 1; i <= dataType.getStateCount(); ++i) {
            if (dataType.getState(Integer.toString(i)) != i - 1) {
                throw new IllegalArgumentException("Data type is not counting numbers.");
            }
        }
    }

    private enum BirthDeathParameterization {
        LINEAR {
            @Override
            double birthRate(int capita, double rate) {
                return capita * rate;
            }

            @Override
            double deathRate(int capita, double rate) {
                return capita * rate;
            }
        };

        abstract double birthRate(int capita, double rate);

        abstract double deathRate(int capita, double rate);
    }
}