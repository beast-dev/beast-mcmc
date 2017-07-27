/*
 * TraitData.java
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

package dr.app.beauti.options;

import dr.stats.DiscreteStatistics;

import java.io.Serializable;
import java.util.*;

/**
 * @author Andrew Rambaut
 */
public class Predictor implements Serializable {
    private static final long serialVersionUID = -9152518508699327745L;

    public enum Type {
        MATRIX,
        ORIGIN_VECTOR,
        DESTINATION_VECTOR,
        BOTH_VECTOR;

        public String toString() {
            return name().toLowerCase();
        }
    }

    private final Type predictorType;
    private final TraitData trait;

    private String name;
    private boolean isIncluded;
    private boolean isLogged;
    private boolean isStandardized;
    private boolean isOrigin;
    private boolean isDestination;

    private final double[][] data;
    private final String[] rowLabels;

    protected final BeautiOptions options;

    public Predictor(BeautiOptions options, String name, TraitData trait, String[] rowLabels, double[][] data, Type predictorType) {
        this.options = options;
        this.name = name;
        this.trait = trait;
        this.rowLabels = rowLabels;
        this.data = data;
        this.predictorType = predictorType;
        this.isIncluded = true;
        this.isLogged = true;
        this.isStandardized = true;
        this.isOrigin = predictorType == Type.ORIGIN_VECTOR || predictorType == Type.BOTH_VECTOR;
        this.isDestination = predictorType == Type.DESTINATION_VECTOR || predictorType == Type.BOTH_VECTOR;
    }

    /////////////////////////////////////////////////////////////////////////

    public Type getType() {
        return predictorType;
    }

    public String getName() {
        return name;
    }

    public void setName(final String name) {
        this.name = name;
    }

    public boolean isIncluded() {
        return isIncluded;
    }

    public void setIncluded(boolean included) {
        isIncluded = included;
    }

    public boolean isLogged() {
        return isLogged;
    }

    public void setLogged(boolean logged) {
        isLogged = logged;
    }

    public boolean isStandardized() {
        return isStandardized;
    }

    public boolean isOrigin() {
        return isOrigin;
    }

    public void setOrigin(boolean origin) {
        isOrigin = origin;
    }

    public boolean isDestination() {
        return isDestination;
    }

    public void setDestination(boolean destination) {
        isDestination = destination;
    }

    public void setStandardized(boolean standardized) {
        isStandardized = standardized;
    }

    public double[][] getMatrixValues(Type predictorType) {
        double[][] matrixValues = new double[data.length][];

        // create a mapping of the states in order the trait has them to that
        // of the imported data.
        List<String> states = new ArrayList<String>(trait.getStatesOfTrait());
        int[] stateIndices = new int[data.length];
        for (int i = 0; i < data.length; i++) {
            stateIndices[states.indexOf(rowLabels[i])] = i;
        }

        if (getType() == Type.MATRIX) {
            if (predictorType != Type.MATRIX) {
                throw new IllegalArgumentException("Predictor is a matrix");
            }
            for (int i = 0; i < data.length; i++) {
                matrixValues[i] = new double[data.length];
                for (int j = 0; j < data.length; j++) {
                    matrixValues[i][j] = data[stateIndices[i]][stateIndices[j]];
                }
            }
        } else {

            if (predictorType != Type.ORIGIN_VECTOR && predictorType != Type.DESTINATION_VECTOR) {
                throw new IllegalArgumentException("Predictor is a vector");
            }

            for (int i = 0; i < data.length; i++) {
                matrixValues[i] = new double[data.length];
            }
            for (int i = 0; i < data.length; i++) {
                for (int j = 0; j < data.length; j++) {
                    if (predictorType == Type.ORIGIN_VECTOR) {
                        matrixValues[i][j] = data[stateIndices[i]][0];
                    } else {
                        matrixValues[j][i] = data[stateIndices[i]][0];
                    }
                }
            }
        }

        // set the diagonal to NaN - this doesn't enter the final vector but helps debugging
        for (int i = 0; i < data.length; i++) {
            matrixValues[i][i] = Double.NaN;
        }

        if (isLogged()) {
            for (int i = 0; i < matrixValues.length; i++) {
                for (int j = 0; j < matrixValues[i].length; j++) {
                    if (i != j) {
                        matrixValues[i][j] = Math.log(matrixValues[i][j]);
                    }
                }
            }
        }

        if (isStandardized()) {
            double[] values = new double[matrixValues.length * (matrixValues.length - 1)];

            int k = 0;
            for (int i = 0; i < matrixValues.length; i++) {
                for (int j = 0; j < matrixValues.length; j++) {

                    if (i != j) {
                        values[k] = matrixValues[i][j];
                        k++;
                    }
                }
            }
            double mean = DiscreteStatistics.mean(values);
            double stdev = DiscreteStatistics.stdev(values);

            for (int i = 0; i < matrixValues.length; i++) {
                for (int j = 0; j < matrixValues.length; j++) {
                    if (i != j) {
                        matrixValues[i][j] = ((matrixValues[i][j] - mean) / stdev);
                    }

                }
            }
        }
        return matrixValues;
    }

    /**
     * Return string with the values in the linear top triangle, bottom triangle format
     * @return
     */
    public String getValueString(Type predictorType) {
        StringBuilder valueString = new StringBuilder();
        double[][] matrix = getMatrixValues(predictorType);

        boolean first = true;

        int n = matrix.length;

        // upper triangle
        for (int i = 0; i < n; i++) {
            for (int j = i + 1; j < n; j++) {
                if (!first) {
                    valueString.append(" ");
                } else {
                    first = false;
                }
                valueString.append(matrix[i][j]);
            }
        }
        // lower triangle
        for (int i = 0; i < n; i++) {
            for (int j = i + 1; j < n; j++) {
                valueString.append(" ");
                valueString.append(matrix[j][i]);
            }
        }
        return valueString.toString();
    }

    public String toString() {
        return name;
    }
}
