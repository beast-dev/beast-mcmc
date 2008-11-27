/*
 * BaseSubstitutionModel.java
 *
 * Copyright (C) 2002-2006 Alexei Drummond and Andrew Rambaut
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

package dr.evomodel.newsubstmodel;

import dr.evolution.datatype.DataType;
import dr.inference.model.AbstractModel;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.ParameterChangeType;

/**
 * An abstract base class for substitution models.
 *
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @version $Id: AbstractSubstitutionModel.java,v 1.41 2005/05/24 20:25:58 rambaut Exp $
 */
@SuppressWarnings({"SuspiciousNameCombination", "UnusedAssignment"})
public abstract class BaseSubstitutionModel extends AbstractModel
        implements SubstitutionModel {


    public static final String MODEL = "model";

    protected DataType dataType = null;

    protected FrequencyModel freqModel;
    protected double[] relativeRates;
    protected double[] storedRelativeRates;

    protected int stateCount;
    protected int rateCount;

    protected boolean eigenInitialised = false;
    protected boolean updateMatrix = true;
    protected boolean storedUpdateMatrix = true;

    private final EigenSystem eigenSystem;

    BaseSubstitutionModel(String name, DataType dataType, FrequencyModel freqModel) {
        this(name, dataType, freqModel, new DefaultEigenSystem(dataType.getStateCount()));
    }

    BaseSubstitutionModel(String name, DataType dataType, FrequencyModel freqModel, EigenSystem eigenSystem) {
        super(name);

        this.eigenSystem = eigenSystem;

        this.dataType = dataType;

        setStateCount(dataType.getStateCount());

        if (freqModel != null) {
            // freqModel can be null at this point but must be
            // in place by the time setupMatrix is called.

            if (freqModel.getDataType() != dataType) {
                throw new IllegalArgumentException("Datatypes do not match!");
            }

            this.freqModel = freqModel;
            addModel(freqModel);

            checkFrequencies();
        }

        updateMatrix = true;
    }

    private void setStateCount(int stateCount) {
        eigenInitialised = false;

        this.stateCount = stateCount;
        rateCount = ((stateCount - 1) * stateCount) / 2;

        relativeRates = new double[rateCount];
        storedRelativeRates = new double[rateCount];
        for (int i = 0; i < rateCount; i++) {
            relativeRates[i] = 1.0;
        }
    }

    // *****************************************************************
    // Interface Model
    // *****************************************************************

    protected void handleModelChangedEvent(Model model, Object object, int index) {
        // frequencyModel changed!
        updateMatrix = true;
        frequenciesChanged();
    }

    protected final void handleParameterChangedEvent(Parameter parameter, int index, ParameterChangeType type) {
        // relativeRates changed
        updateMatrix = true;
        ratesChanged();
    }

    protected void storeState() {

        storedUpdateMatrix = updateMatrix;

        System.arraycopy(relativeRates, 0, storedRelativeRates, 0, rateCount);

        System.arraycopy(Eval, 0, storedEval, 0, stateCount);
        for (int i = 0; i < stateCount; i++) {
            System.arraycopy(Ievc[i], 0, storedIevc[i], 0, stateCount);
            System.arraycopy(Evec[i], 0, storedEvec[i], 0, stateCount);
        }

    }

    /**
     * Restore the additional stored state
     */
    protected void restoreState() {

        updateMatrix = storedUpdateMatrix;

        // To restore all this stuff just swap the pointers...
        double[] tmp1 = storedRelativeRates;
        storedRelativeRates = relativeRates;
        relativeRates = tmp1;

        tmp1 = storedEval;
        storedEval = Eval;
        Eval = tmp1;

        double[][] tmp2 = storedIevc;
        storedIevc = Ievc;
        Ievc = tmp2;

        tmp2 = storedEvec;
        storedEvec = Evec;
        Evec = tmp2;

    }

    protected void acceptState() {
    } // nothing to do

    abstract protected void frequenciesChanged();

    abstract protected void ratesChanged();

    abstract protected void setupRelativeRates(double[] rates);

    public FrequencyModel getFrequencyModel() {
        return freqModel;
    }

    /**
     * @return the data type
     */
    public DataType getDataType() {
        return dataType;
    }

    /**
     * get the complete transition probability matrix for the given distance
     *
     * @param distance the expected number of substitutions
     * @param matrix   an array to store the matrix
     */
    public void getTransitionProbabilities(double distance, double[] matrix) {
        double temp;

        // this must be synchronized to avoid being called simultaneously by
        // two different likelihood threads - AJD
        synchronized (this) {
            if (updateMatrix) {
                setupMatrix();
            }
        }

        double[][] Evec = eigenSystem.getEigenVectors();
//        if (DEBUG) System.err.println(new dr.math.matrixAlgebra.Vector(Evec[0]));
        double[][] Ievc = eigenSystem.getInverseEigenVectors();
//        if (DEBUG) System.err.println(new dr.math.matrixAlgebra.Vector(Ievc[0]));

        // implemented a pool of iexp matrices to support multiple threads
        // without creating a new matrix each call. - AJD
        double[][] iexp = new double[stateCount][stateCount];
        for (int i = 0; i < stateCount; i++) {
            temp = Math.exp(distance * Eval[i]);
            for (int j = 0; j < stateCount; j++) {
                iexp[i][j] = Ievc[i][j] * temp;
            }
        }

        int u = 0;
        for (int i = 0; i < stateCount; i++) {
            for (int j = 0; j < stateCount; j++) {
                temp = 0.0;
                for (int k = 0; k < stateCount; k++) {
                    temp += Evec[i][k] * iexp[k][j];
                }

                matrix[u] = Math.abs(temp);
                u++;
            }
        }
    }

    /**
     * This function returns the Eigen vectors.
     *
     * @return the array
     */
    public double[][] getEigenVectors() {
        synchronized (this) {
            if (updateMatrix) {
                setupMatrix();
            }
        }
        return Evec;
    }

    /**
     * This function returns the inverse Eigen vectors.
     *
     * @return the array
     */
    public double[][] getInverseEigenVectors() {
        synchronized (this) {
            if (updateMatrix) {
                setupMatrix();
            }
        }
        return Ievc;
    }

    /**
     * This function returns the Eigen values.
     */
    public double[] getEigenValues() {
        synchronized (this) {
            if (updateMatrix) {
                setupMatrix();
            }
        }
        return Eval;
    }

    /**
     * setup substitution matrix
     */
    protected void setupMatrix() {
        setupRelativeRates(relativeRates);

        int k = 0;
        // Set the instantaneous rate matrix
        for (int i = 0; i < stateCount; i++) {
            for (int j = i + 1; j < stateCount; j++) {
                q[i][j] = relativeRates[k] * freqModel.getFrequency(j);
                q[j][i] = relativeRates[k] * freqModel.getFrequency(i);
                k += 1;
            }
        }
        makeValid(q, stateCount);
        normalize(q, freqModel.getFrequencies());

        updateMatrix = false;
    }

    // Make it a valid rate matrix (make sum of rows = 0)
    void makeValid(double[][] matrix, int dimension) {
        for (int i = 0; i < dimension; i++) {
            double sum = 0.0;
            for (int j = 0; j < dimension; j++) {
                if (i != j)
                    sum += matrix[i][j];
            }
            matrix[i][i] = -sum;
        }
    }

    /**
     * Normalize rate matrix to one expected substitution per unit time
     *
     * @param matrix the matrix to normalize to one expected substitution
     * @param pi     the equilibrium distribution of states
     */
    void normalize(double[][] matrix, double[] pi) {
        double subst = 0.0;
        int dimension = pi.length;

        for (int i = 0; i < dimension; i++)
            subst += -matrix[i][i] * pi[i];

        for (int i = 0; i < dimension; i++) {
            for (int j = 0; j < dimension; j++) {
                matrix[i][j] = matrix[i][j] / subst;
            }
        }
    }

    /**
     * Ensures that frequencies are not smaller than MINFREQ and
     * that two frequencies differ by at least 2*MINFDIFF.
     * This avoids potential problems later when eigenvalues
     * are computed.
     */
    private void checkFrequencies() {
        // required frequency difference
        double MINFDIFF = 1.0E-10;

        // lower limit on frequency
        double MINFREQ = 1.0E-10;

        int maxi = 0;
        double sum = 0.0;
        double maxfreq = 0.0;
        for (int i = 0; i < stateCount; i++) {
            double freq = freqModel.getFrequency(i);
            if (freq < MINFREQ) freqModel.setFrequency(i, MINFREQ);
            if (freq > maxfreq) {
                maxfreq = freq;
                maxi = i;
            }
            sum += freqModel.getFrequency(i);
        }
        double diff = 1.0 - sum;
        freqModel.setFrequency(maxi, freqModel.getFrequency(maxi) + diff);

        for (int i = 0; i < stateCount - 1; i++) {
            for (int j = i + 1; j < stateCount; j++) {
                if (freqModel.getFrequency(i) == freqModel.getFrequency(j)) {
                    freqModel.setFrequency(i, freqModel.getFrequency(i) + MINFDIFF);
                    freqModel.setFrequency(j, freqModel.getFrequency(j) - MINFDIFF);
                }
            }
        }
    }

    /**
     * allocate memory for the Eigen routines
     */
    protected void initialiseEigen() {

        Eval = new double[stateCount];
        Evec = new double[stateCount][stateCount];
        Ievc = new double[stateCount][stateCount];

        storedEval = new double[stateCount];
        storedEvec = new double[stateCount][stateCount];
        storedIevc = new double[stateCount][stateCount];

        q = new double[stateCount][stateCount];

    }

    double q[][];

    // Eigenvalues, eigenvectors, and inverse eigenvectors
    protected double[] Eval;
    protected double[] storedEval;
    protected double[][] Evec;
    protected double[][] storedEvec;
    protected double[][] Ievc;
    protected double[][] storedIevc;


}