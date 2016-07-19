/*
 * TwoStateCovarionModel.java
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

package dr.oldevomodel.substmodel;

import dr.evolution.datatype.TwoStateCovarion;
import dr.oldevomodelxml.substmodel.TwoStateCovarionModelParser;
import dr.inference.model.Parameter;

/**
 * @author Helen Shearman
 * @author Alexei Drummond
 * @version $Id$
 */
public class TwoStateCovarionModel extends AbstractSubstitutionModel {

    /**
     * constructor
     *
     * @param dataType           the data type
     * @param freqModel          the frequencies
     * @param alphaParameter     the rate of evolution in slow mode
     * @param switchingParameter the rate of flipping between slow and fast modes
     */
    public TwoStateCovarionModel(TwoStateCovarion dataType, FrequencyModel freqModel,
                                 Parameter alphaParameter,
                                 Parameter switchingParameter) {
        super(TwoStateCovarionModelParser.COVARION_MODEL, dataType, freqModel);

        alpha = alphaParameter;
        this.switchingParameter = switchingParameter;

        addVariable(alpha);
        addVariable(switchingParameter);
        setupRelativeRates();
    }

    protected void frequenciesChanged() {
        // DO NOTHING
    }

    protected void ratesChanged() {
        setupRelativeRates();
    }

    protected void setupRelativeRates() {

        relativeRates[0] = alpha.getParameterValue(0);
        relativeRates[1] = switchingParameter.getParameterValue(0);
        relativeRates[2] = 0.0;
        relativeRates[3] = 0.0;
        relativeRates[4] = switchingParameter.getParameterValue(0);
        relativeRates[5] = 1.0;
    }

    public String toString() {

        return SubstitutionModelUtils.toString(relativeRates, dataType, true, 2);
    }

    /**
     * Normalize rate matrix to one expected substitution per unit time
     *
     * @param matrix the matrix to normalize to one expected substitution
     * @param pi     the equilibrium distribution of states
     */
    void normalize(double[][] matrix, double[] pi) {

        if (isNormalized) {

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

            //System.out.println("switchingProportion=" + switchingProportion);

            // normalize, removing switches
            for (int i = 0; i < dimension; i++) {
                for (int j = 0; j < dimension; j++) {
                    matrix[i][j] = matrix[i][j] / (1.0 - switchingProportion);
                }
            }
        }
    }

    private Parameter alpha;
    private Parameter switchingParameter;

    // if true then matrix will be normalized to output 1 substitution per unit time
    private boolean isNormalized = true;

}
