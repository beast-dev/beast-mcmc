/*
 * BinaryCovarionModel.java
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

import dr.evomodelxml.substmodel.BinaryCovarionModelParser;
import dr.evolution.datatype.TwoStateCovarion;
import dr.oldevomodel.substmodel.SubstitutionModelUtils;
import dr.inference.model.Parameter;

/**
 * @author Alexei Drummond
 * @author Marc A. Suchard
 *         <p/>
 *         a	the rate of the slow rate class
 *         1	the rate of the fast rate class
 *         p0	the equilibrium frequency of zero states
 *         p1	1 - p0, the equilibrium frequency of one states
 *         f0	the equilibrium frequency of slow rate class
 *         f1	1 - f0, the equilibrium frequency of fast rate class
 *         s	the rate of switching
 *         <p/>
 *         then the (unnormalized) instantaneous rate matrix (unnormalized Q) should be
 *         <p/>
 *         [ -(a*p1)-s ,   a*p1    ,    s   ,   0   ]
 *         [   a*p0    , -(a*p0)-s ,    0   ,   s   ]
 *         [    s      ,     0     ,  -p1-s ,  p1   ]
 *         [    0      ,     s     ,    p0  , -p0-s ]
 */
public class BinaryCovarionModel extends AbstractCovarionModel {

    /**
     * @param dataType           the data type
     * @param frequencies        the frequencies of the visible states
     * @param hiddenFrequencies  the frequencies of the hidden rates
     * @param alphaParameter     the rate of evolution in slow mode
     * @param switchingParameter the rate of flipping between slow and fast modes
     */
    public BinaryCovarionModel(TwoStateCovarion dataType,
                               Parameter frequencies,
                               Parameter hiddenFrequencies,
                               Parameter alphaParameter,
                               Parameter switchingParameter,
                               dr.oldevomodel.substmodel.BinaryCovarionModel.Version version) {

        super(BinaryCovarionModelParser.COVARION_MODEL, dataType, frequencies, hiddenFrequencies);

        this.alpha = alphaParameter;
        this.switchRate = switchingParameter;
        this.frequencies = frequencies;
        this.hiddenFrequencies = hiddenFrequencies;
        this.version = version;

        addVariable(alpha);
        addVariable(switchRate);
        addVariable(frequencies);
        addVariable(hiddenFrequencies);
    }

    protected void setupQMatrix(double[] rates, double[] pi, double[][] unnormalizedQ) {

        double a = alpha.getParameterValue(0);
        double s = switchRate.getParameterValue(0);
        double f0 = hiddenFrequencies.getParameterValue(0);
        double f1 = hiddenFrequencies.getParameterValue(1);
        double p0 = frequencies.getParameterValue(0);   // TODO Use pi[0], ..., pi[3]
        double p1 = frequencies.getParameterValue(1);

        assert Math.abs(1.0 - f0 - f1) < 1e-8;
        assert Math.abs(1.0 - p0 - p1) < 1e-8;

        f0 = version.getF0(f0);
        f1 = version.getF1(f1);

        unnormalizedQ[0][1] = a * p1;
        unnormalizedQ[0][2] = s * f0;
        unnormalizedQ[0][3] = 0.0;

        unnormalizedQ[1][0] = a * p0;
        unnormalizedQ[1][2] = 0.0;
        unnormalizedQ[1][3] = s * f0;

        unnormalizedQ[2][0] = s * f1;
        unnormalizedQ[2][1] = 0.0;
        unnormalizedQ[2][3] = p1;

        unnormalizedQ[3][0] = 0.0;
        unnormalizedQ[3][1] = s * f1;
        unnormalizedQ[3][2] = p0;
    }

    protected void frequenciesChanged() {
    }

    protected void ratesChanged() {
    }

    public String toString() {

        final int stateCount = dataType.getStateCount();
        double[] tmp = new double[stateCount * stateCount];
        getInfinitesimalMatrix(tmp);
        double[][] q = new double[stateCount][stateCount];
        for (int i = 0; i < stateCount; ++i) {
            for (int j = 0; j < stateCount; ++j) {
                q[i][j] = tmp[i * stateCount + j];
            }
        }
        setupQMatrix(null, null, q);
        return SubstitutionModelUtils.toString(q, dataType, 2);
    }

    /**
     * Normalize rate matrix to one expected substitution per unit time
     *
     * @param matrix the matrix to normalize to one expected substitution
     * @param pi     the equilibrium distribution of states
     */
    protected double getNormalizationValue(double[][] matrix, double[] pi) {

        double subst = 0.0;
        int dimension = pi.length;

        for (int i = 0; i < dimension; i++) {
            subst += -matrix[i][i] * pi[i];
        }

        // normalize, including switches
        for (int i = 0; i < dimension; i++) {
            for (int j = 0; j < dimension; j++) {
                matrix[i][j] = matrix[i][j] / subst;
            }
        }

        double switchingProportion = 0.0;
        switchingProportion += matrix[0][2] * pi[2];
        switchingProportion += matrix[2][0] * pi[0];
        switchingProportion += matrix[1][3] * pi[3];
        switchingProportion += matrix[3][1] * pi[1];

        // normalize, removing switches
        for (int i = 0; i < dimension; i++) {
            for (int j = 0; j < dimension; j++) {
                matrix[i][j] = matrix[i][j] / (1.0 - switchingProportion);
            }
        }
        return 1.0; // Already normalized
    }

    private Parameter alpha;
    private Parameter switchRate;
    private Parameter frequencies;
    private Parameter hiddenFrequencies;
    final private dr.oldevomodel.substmodel.BinaryCovarionModel.Version version;

}