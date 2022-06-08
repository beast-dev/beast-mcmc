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

import java.util.List;

/**
 * @author Marc Suchard
 */
public class BirthDeathSubstitutionModel extends ComplexSubstitutionModelAtStationarity {

    private final List<Parameter> parameters;
    private final BirthDeathParameterization parameterization;
    private final boolean useStationaryDistribution;
    private final int dim;

    public BirthDeathSubstitutionModel(String name,
                                       List<Parameter> parameters,
                                       DataType dataType,
                                       boolean useStationaryDistribution) {
        super(name, dataType, null);

        this.parameters = parameters;

        int len = 0;
        for (Parameter p : parameters) {
            addVariable(p);
            len += p.getDimension();
        }
        this.dim = len;

        this.useStationaryDistribution = useStationaryDistribution;
        this.parameterization = BirthDeathParameterization.LINEAR;
        this.freqModel = setupEquilibriumModel();

        checkDataType(dataType);

        if (dim != parameterization.getRequiredDimension()) {
            throw new IllegalArgumentException("Invalid parameterization for birth-death substitution process");
        }
    }

    @Override
    protected void frequenciesChanged() {
        throw new RuntimeException("Frequencies are fixed to birth-death equilibrium distribution");
    }

    @Override
    protected double[] getPi() {
        return freqModel.getFrequencies();
    }

    @Override
    protected void setupRelativeRates(double[] rates) {
        // Do nothing
    }

    @Override
    protected void setupQMatrix(double[] rates, double[] pi, double[][] matrix) {

        final double[] perCapitaRate = getRates();

        for (int i = 1; i < stateCount; i++) {
            matrix[i - 1][i] = parameterization.birthRate(i,  perCapitaRate);
            matrix[i][i - 1] = parameterization.deathRate(i + 1,  perCapitaRate);
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

    private double[] getRates() {

        double[] result = new double[dim];

        int offset = 0;
        for (Parameter p : parameters) {
            for (int i = 0; i < p.getDimension(); ++i) {
                result[offset + i] = p.getParameterValue(i);
            }
            offset += p.getDimension();
        }

        return result;
    }
    
    private FrequencyModel setupEquilibriumModel() {
        Parameter equilibrium = new Parameter.Proxy("equilibrium", stateCount) {

            @Override
            public double getParameterValue(int dim) {
                if (useStationaryDistribution) {
                    return getStationaryDistribution()[dim];
                } else {
                    return 1.0 / stateCount;
                }
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

    // Some different models:
    //   https://doi.org/10.1093/molbev/msg084
    //   https://doi.org/10.1534/genetics.103.022665
    //   https://doi.org/10.1534/genetics.110.125260

    private enum BirthDeathParameterization {
        LINEAR {
            @Override
            double birthRate(int capita, double[] rate) {
                return capita * rate[0];
            }

            @Override
            double deathRate(int capita, double[] rate) {
                return capita * rate[1];
            }

            @Override
            int getRequiredDimension() {
                return 2;
            }
        },
        QUADRATIC {
            @Override
            double birthRate(int capita, double[] rate) {
                return capita * rate[0]
                        + capita * capita * rate[1];
            }

            @Override
            double deathRate(int capita, double[] rate) {
                return capita * rate[2]
                        + capita * capita * rate[3];
            }

            @Override
            int getRequiredDimension() {
                return 4;
            }
        };

        abstract double birthRate(int capita, double[] rate);

        abstract double deathRate(int capita, double[] rate);

        abstract int getRequiredDimension();
    }
}