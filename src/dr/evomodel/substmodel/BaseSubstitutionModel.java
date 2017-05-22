/*
 * BaseSubstitutionModel.java
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

package dr.evomodel.substmodel;

import dr.evolution.datatype.DataType;
import dr.inference.model.AbstractModel;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;

import java.util.Arrays;

/**
 * An abstract base class for substitution models.
 *
 * @author Andrew Rambaut
 * @author Marc Suchard
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

    protected boolean updateMatrix = true;
    protected boolean storedUpdateMatrix = true;

    private final EigenSystem eigenSystem;

    public BaseSubstitutionModel(String name) {
        super(name);

        // For a wrapper model (KroneckerSumSM), most computation is handled in the wrapped classes
        eigenSystem = null;
        q = null;
    }

    public BaseSubstitutionModel(String name, DataType dataType, FrequencyModel freqModel) {
        this(name, dataType, freqModel, null);
    }

    public BaseSubstitutionModel(String name, DataType dataType, FrequencyModel freqModel, EigenSystem eigenSystem) {
        super(name);

        if (eigenSystem == null)
            this.eigenSystem = getDefaultEigenSystem(dataType.getStateCount());
        else
            this.eigenSystem = eigenSystem;

        this.dataType = dataType;

        setStateCount(dataType.getStateCount());

        if (freqModel != null) {
            // freqModel can be null at this point but must be
            // in place by the time setupMatrix is called.

            if (freqModel.getDataType() != dataType) {
                throw new IllegalArgumentException("Datatypes do not match.");
            }

            this.freqModel = freqModel;
            addModel(freqModel);

            if (!(freqModel instanceof CovarionFrequencyModel)) {
                checkFrequencies();
            }
        }

        q = new double[stateCount][stateCount];

        updateMatrix = true;
    }

    protected EigenSystem getDefaultEigenSystem(int stateCount) {
        return new DefaultEigenSystem(stateCount);
    }

    private void setStateCount(int stateCount) {
        this.stateCount = stateCount;
        rateCount = getRateCount(stateCount);

        relativeRates = new double[rateCount];
        storedRelativeRates = new double[rateCount];
        for (int i = 0; i < rateCount; i++) {
            relativeRates[i] = 1.0;
        }
    }

    protected int getRateCount(int stateCount) {
        return ((stateCount - 1) * stateCount) / 2;
    }

    // *****************************************************************
    // Interface Model
    // *****************************************************************

    protected void handleModelChangedEvent(Model model, Object object, int index) {
        // frequencyModel changed!
        updateMatrix = true;
        frequenciesChanged();
        fireModelChanged();
    }

    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        // relativeRates changed
        updateMatrix = true;
        ratesChanged();
    }

    protected void storeState() {

        storedUpdateMatrix = updateMatrix;

        System.arraycopy(relativeRates, 0, storedRelativeRates, 0, rateCount);

        if (eigenDecomposition != null) {
            storedEigenDecomposition = eigenDecomposition.copy();
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

        EigenDecomposition tmp = storedEigenDecomposition;
        storedEigenDecomposition = eigenDecomposition;
        eigenDecomposition = tmp;

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

        EigenDecomposition eigen = getEigenDecomposition();

        if (eigen == null) {
            Arrays.fill(matrix, 0.0);
            return;
        }

        double[] Evec = eigen.getEigenVectors();
        double[] Ievc = eigen.getInverseEigenVectors();
        double[] Eval = eigen.getEigenValues();

        // implemented a pool of iexp matrices to support multiple threads
        // without creating a new matrix each call. - AJD
        double[][] iexp = new double[stateCount][stateCount];
        for (int i = 0; i < stateCount; i++) {
            temp = Math.exp(distance * Eval[i]);
            for (int j = 0; j < stateCount; j++) {
                iexp[i][j] = Ievc[i * stateCount + j] * temp;
            }
        }

        int u = 0;
        for (int i = 0; i < stateCount; i++) {
            for (int j = 0; j < stateCount; j++) {
                temp = 0.0;
                for (int k = 0; k < stateCount; k++) {
                    temp += Evec[i * stateCount + k] * iexp[k][j];
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
    public EigenDecomposition getEigenDecomposition() {
        synchronized (this) {
            if (updateMatrix) {
                decompose();
            }
        }
        return eigenDecomposition;
    }

    protected void setupQMatrix(double[] rates, double[] pi, double[][] matrix) {
        int k = 0;
        // Set the instantaneous rate matrix
        for (int i = 0; i < stateCount; i++) {
            for (int j = i + 1; j < stateCount; j++) {
                matrix[i][j] = rates[k] * pi[j];
                matrix[j][i] = rates[k] * pi[i];
                k++;
            }
        }
    }

    /**
     * setup substitution matrix
     */
    private void decompose() {

        double normalization = setupMatrix();

        eigenDecomposition = eigenSystem.decomposeMatrix(q);

        if (eigenDecomposition != null)
            eigenDecomposition.normalizeEigenValues(normalization);

        updateMatrix = false;
    }

    private double setupMatrix() {
        setupRelativeRates(relativeRates);
        double[] pi = freqModel.getFrequencies();
        setupQMatrix(relativeRates, pi, q);
        makeValid(q, stateCount);
        return getNormalizationValue(q, pi);
    }

    public void getInfinitesimalMatrix(double[] out) {

        double normalization = setupMatrix();
        int index = 0;
        for (int i = 0; i < stateCount; i++) {
            for (int j = 0; j < stateCount; j++) {
                out[index] = q[i][j] / normalization;
                index++;
            }
        }
    }

    /**
     * Normalize rate matrix to one expected substitution per unit time
     *
     * @param matrix the matrix to normalize to one expected substitution
     * @param pi     the equilibrium distribution of states
     * @return expected substitution per unit time
     */
    protected double getNormalizationValue(double[][] matrix, double[] pi) {
        double subst = 0.0;

        for (int i = 0; i < stateCount; i++)
            subst += -matrix[i][i] * pi[i];

        return subst;
    }

    static String format = "%2.1e";

    public String printQ() {
        double[] out = new double[stateCount * stateCount];
        getInfinitesimalMatrix(out);

        StringBuffer sb = new StringBuffer();
        int index = 0;
        for (int i = 0; i < stateCount; i++) {
            if (i > 0)
                sb.append(String.format(format, out[index]));
            for (int j = 1; j < stateCount; j++) {
                sb.append("\t");
                if (j != i)
                    sb.append(String.format(format, out[index]));
                index++;
            }
            sb.append("\n");
        }
        return sb.toString();
    }

    // Make it a valid rate matrix (make sum of rows = 0)
    private void makeValid(double[][] matrix, int dimension) {
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
     * Ensures that frequencies are not smaller than MINFREQ and
     * that two frequencies differ by at least 2*MINFDIFF.
     * This avoids potential problems later when eigenvalues
     * are computed.
     */
    private void checkFrequencies() {
        // required frequency difference
        double MINFDIFF = getMINFDIFF();

        // lower limit on frequency
        double MINFREQ = getMINFREQ();


//        System.err.println("diff = "+MINFDIFF+" freq = "+MINFREQ);
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

    public boolean canReturnComplexDiagonalization() {
        return false;
    }

    protected double getMINFDIFF() {
        return 1.0E-10;
    }

    protected double getMINFREQ() {
        return 1.0E-10;
    }

    protected double[][] getQCopy() {
        double[][] copy = new double[q.length][q.length];
        for (int i = 0; i < q.length; ++i) {
            System.arraycopy(q[i], 0, copy[i], 0, q.length);
        }
        return copy;
    }

    private final double q[][];
    protected EigenDecomposition eigenDecomposition;
    private EigenDecomposition storedEigenDecomposition;

}