/*
 * SVSComplexSubstitutionModel.java
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

import dr.evolution.datatype.DataType;
import dr.inference.model.BayesianStochasticSearchVariableSelection;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;
import dr.math.matrixAlgebra.Vector;

/**
 * @author Marc A. Suchard
 */
public class SVSComplexSubstitutionModel extends ComplexSubstitutionModel implements BayesianStochasticSearchVariableSelection {

    public SVSComplexSubstitutionModel(String name, DataType dataType, FrequencyModel rootFreqModel,
                                       Parameter rates, Parameter indicators) {
        super(name, dataType, rootFreqModel, rates);
        if (indicators != null) {
            this.indicators = indicators;
            addVariable(indicators);
        } else {
            this.indicators = new Parameter.Default(rates.getDimension(), 1.0);
        }
    }

    protected double[] getRates() {
        double[] rates = infinitesimalRates.getParameterValues();
        double[] indies = indicators.getParameterValues();
        final int dim = rates.length;
        for (int i = 0; i < dim; i++)
            rates[i] *= indies[i];
        return rates;
    }

    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        if (variable == infinitesimalRates && indicators.getParameterValue(index) == 0)
            return;  // ignore, does not affect likelihood
        super.handleVariableChangedEvent(variable, index, type);
    }

    public Parameter getIndicators() {
        return indicators;
    }

//     public double getLogLikelihood() {
//        double logL = super.getLogLikelihood();
//        if (logL == 0 && considerConnectedness && !isStronglyConnected()) { // Also check that graph is connected
//            logL = Double.NEGATIVE_INFINITY;
//        }
//        return logL;
//    }

//    private boolean hasEdge(int i, int j) {
//        return i != j && getIndicators().getParameterValue(getEntry(i,j,stateCount)) == 1;
//    }
//
//    private static int getEntry(int i, int j, int dim) {
//        int entry = i * (dim - 1) + j;
//        if( j > i )
//            entry--;
//        return entry;
//    }
//
//    private void depthFirstSearch(int node, BitVector visited) {
//        visited.set(node);
//        for(int v=0; v<stateCount; v++) {
//            if (hasEdge(node,v) && !visited.get(v))
//                depthFirstSearch(v,visited);
//        }
//    }

//    /* Determines if the graph is strongly connected, such that there exists
//     * a directed path from any vertex to any other vertex
//     *
//     */
//    public boolean isStronglyConnected() {
//        BitVector visited = new BitVector(stateCount);
//        boolean connected = true;
//        for(int i=0; i<stateCount && connected; i++) {
//            visited.clear();
//            depthFirstSearch(i,visited);
//            connected = visited.cardinality() == stateCount;
//        }
//        return connected;
//    }

    public boolean validState() {
        return getLogLikelihood() == 0;
    }

    public static void main(String[] arg) {

        double[] rates = new double[]{0.097,0.515,3.346,0.623,0.389,0.631,0.362,1.127,4.262,0.424,0.758,1.297,0.728,0.228,0.003,0.075,0.312,0.356,0.927,4.420,2.719,0.264,4.267,0.741,1.106,1.568,1.215,0.172,0.204,1.493,0.592,0.105,1.583,1.201,0.783,2.224,0.888,1.401,0.137,0.259,1.197,0.056,1.939,1.385,0.690,0.815,0.279,1.100,1.715,0.011,1.509,0.961,0.112,1.305,2.797,0.578,1.177,1.009,0.316,1.143,1.861,0.176,0.140,0.104,0.571,0.521,0.761,1.795,1.065,1.563,2.972,2.295,0.075,1.690,1.011,0.128,0.484,0.355,1.668,1.052,0.089,0.104,0.056,1.591,0.054,0.487,1.034,1.145,0.403,0.254,0.474,0.181,0.124,2.067,2.208,0.120,2.638,0.195,1.897,1.955,1.113,1.399,4.901,3.218,0.361,1.934,1.681,3.572,0.806,0.077,0.042,0.310,0.243,1.526,3.572,4.173,0.740,3.086,0.645,2.755,0.280,2.476,0.476,5.174,0.057,0.225,1.310,0.201,0.491,1.507,0.604,0.404,0.907,0.048,2.439,0.676,0.382,0.062,1.173,2.026,2.813,0.655,1.511,0.452,0.137,0.435,0.685,0.826,1.735,0.935,0.697,1.230,5.416,1.043,0.747,0.945};
        double[] indicators = new double[]{1.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,1.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,1.000,0.000,0.000,0.000,0.000,1.000,0.000,0.000,0.000,0.000,1.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,1.000,0.000,0.000,0.000,0.000,0.000,1.000,0.000,0.000,1.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,1.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,1.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,1.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,1.000,0.000,0.000,0.000,0.000,1.000,0.000,0.000,0.000,0.000,0.000,1.000,0.000,0.000,0.000,0.000,1.000,0.000,0.000,0.000,1.000,1.000,1.000,0.000,1.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,1.000};

        Parameter ratesP = new Parameter.Default(rates);
        Parameter indicatorsP = new Parameter.Default(indicators);

        DataType dataType = new DataType() {

            public String getDescription() {
                return null;
            }

            public int getType() {
                return 0;
            }

            @Override
            public char[] getValidChars() {
                return null;
            }

            public int getStateCount() {
                return 13;
            }
        };

        FrequencyModel freqModel = new FrequencyModel(dataType, new Parameter.Default(400, 1.0/400.0));

        SVSComplexSubstitutionModel substModel = new SVSComplexSubstitutionModel("test",
                  dataType,
                  freqModel,
                  ratesP,indicatorsP);
        double logL = substModel.getLogLikelihood();
        System.out.println("Prior = "+logL);
        if( !Double.isInfinite(logL) ) {
            double[] finiteTimeProbs = new double[substModel.getDataType().getStateCount() * substModel.getDataType().getStateCount()];
            double time = 1.0;
            substModel.getTransitionProbabilities(time,finiteTimeProbs);
            System.out.println("Probs = "+new Vector(finiteTimeProbs));
        }
    }
    
    private Parameter indicators;
//    private boolean considerConnectedness = false;

}
