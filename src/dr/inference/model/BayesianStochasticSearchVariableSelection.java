/*
 * BayesianStochasticSearchVariableSelection.java
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

package dr.inference.model;

import cern.colt.bitvector.BitVector;
import dr.math.MathUtils;
import dr.oldevomodel.substmodel.SubstitutionModel;

/**
 * @author Marc Suchard
 */

public interface BayesianStochasticSearchVariableSelection {

    public Parameter getIndicators();

    public boolean validState();

    public class Utils {

        public static boolean connectedAndWellConditioned(double[] probability,
                                                          dr.evomodel.substmodel.SubstitutionModel substModel) {
            if (probability == null) {
                int stateCount = substModel.getDataType().getStateCount();
                probability = new double[stateCount*stateCount];
            }
            try {
                substModel.getTransitionProbabilities(defaultExpectedMutations,probability);
                return connectedAndWellConditioned(probability);
            } catch (Exception e) { // Any numerical error is bad news
                return false;
            }
        }

        public static boolean connectedAndWellConditioned(double[] probability,
                                                          SubstitutionModel substModel) {
            if (probability == null) {
                int stateCount = substModel.getDataType().getStateCount();
                probability = new double[stateCount*stateCount];
            }
            try {
                substModel.getTransitionProbabilities(defaultExpectedMutations,probability);
                return connectedAndWellConditioned(probability);
            } catch (Exception e) { // Any numerical error is bad news
                return false;
            }
        }
                
        public static boolean connectedAndWellConditioned(double[] probability) {
            for(double prob : probability) {
                if(prob < tolerance || prob >= 1.0) {                    
                    return false;
                }
            }
            return true;
        }

        public static void randomize(Parameter indicators,int dim, boolean reversible) {
            do {
                for (int i = 0; i < indicators.getDimension(); i++)
                    indicators.setParameterValue(i,
                            (MathUtils.nextDouble() < 0.5) ? 0.0 : 1.0);
            } while (!(isStronglyConnected(indicators.getParameterValues(),
                    dim, reversible)));
        }

        public static void setTolerance(double newTolerance) {
            tolerance = newTolerance;
        }

        public static double getTolerance() {
            return tolerance;
        }

        public static void setScalar(double newScalar) {
            defaultExpectedMutations = newScalar;
        }
        public static double getScalar() {
            return defaultExpectedMutations;
        }

        /* Determines if the graph is strongly connected, such that there exists
        * a directed path from any vertex to any other vertex
        *
        */
        public static boolean isStronglyConnected(double[] indicatorValues, int dim, boolean reversible) {
            BitVector visited = new BitVector(dim);
            boolean connected = true;
            for (int i = 0; i < dim && connected; i++) {
                visited.clear();
                depthFirstSearch(i, visited, indicatorValues, dim, reversible);
                connected = visited.cardinality() == dim;
            }
            return connected;
        }

        private static boolean hasEdge(int i, int j, double[] indicatorValues,
                                      int dim, boolean reversible) {
            return i != j && indicatorValues[getEntry(i, j, dim, reversible)] == 1;
        }

        private static int getEntry(int i, int j, int dim, boolean reversible) {
            if (reversible) {
                if (j < i) {
                    return getEntry(j,i,dim,reversible);
                }
                int entry = i * dim - i * (i + 1) / 2 + j - 1 -i;
                return entry;
            }

            int entry = i * (dim - 1) + j;
            if (j > i)
                entry--;
            return entry;
        }

        private static void depthFirstSearch(int node, BitVector visited, double[] indicatorValues,
                                            int dim, boolean reversible) {
            visited.set(node);
            for (int v = 0; v < dim; v++) {
                if (hasEdge(node, v, indicatorValues, dim, reversible) && !visited.get(v))
                    depthFirstSearch(v, visited, indicatorValues, dim, reversible);
            }
        }

        private static double defaultExpectedMutations = 1.0;
        private static double tolerance = 1E-20;
    }
}
