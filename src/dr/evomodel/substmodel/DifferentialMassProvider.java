/*
 * DifferentialMassProvider.java
 *
 * Copyright (c) 2002-2018 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

import dr.math.matrixAlgebra.WrappedMatrix;

/**
 * @author Marc A. Suchard
 * @author Xiang Ji
 */
public interface DifferentialMassProvider {

    double[] getDifferentialMassMatrix(double time);

    enum Mode {
        EXACT {
            @Override
            public double[] dispatch(double time,
                              DifferentiableSubstitutionModel model,
                              WrappedMatrix infinitesimalDifferentialMatrix) {
               
                return DifferentiableSubstitutionModelUtil.getExactDifferentialMassMatrix(
                        time, model.getDataType().getStateCount(),
                        infinitesimalDifferentialMatrix, model.getEigenDecomposition());
            }

            @Override
            public String getReport() {
                return "Exact";
            }
        },
        APPROXIMATE {
            @Override
            public double[] dispatch(double time,
                              DifferentiableSubstitutionModel model,
                              WrappedMatrix infinitesimalDifferentialMatrix) {

                return DifferentiableSubstitutionModelUtil.getApproximateDifferentialMassMatrix(
                        time, model.getDataType().getStateCount(),
                        infinitesimalDifferentialMatrix, model.getEigenDecomposition());
            }

            @Override
            public String getReport() {
                return "Approximate wrt parameter";
            }
        },
        AFFINE {
            @Override
            public double[] dispatch(double time,
                              DifferentiableSubstitutionModel model,
                              WrappedMatrix infinitesimalDifferentialMatrix) {

                throw new RuntimeException("Not yet implemented");
            }

            @Override
            public String getReport() {
                return "Affine-corrected wrt parameter";
            }
        };

        public abstract double[] dispatch(double time,
                                   DifferentiableSubstitutionModel model,
                                   WrappedMatrix infinitesimalDifferentialMatrix);

        public abstract String getReport();
    }

    class DifferentialWrapper implements DifferentialMassProvider {
        
        private final DifferentiableSubstitutionModel baseModel;
        private final WrtParameter wrt;
        private final Mode mode;

        public DifferentialWrapper(DifferentiableSubstitutionModel baseModel,
                                   WrtParameter wrt,
                                   Mode mode) {
            this.baseModel = baseModel;
            this.wrt = wrt;
            this.mode = mode;
        }

        @Override
        public double[] getDifferentialMassMatrix(double time) {
            return mode.dispatch(time, baseModel, baseModel.getInfinitesimalDifferentialMatrix(wrt));
        }

        public interface WrtParameter {

            double getRate(int switchCase);

            double getNormalizationDifferential();

            void setupDifferentialFrequencies(double[] differentialFrequencies, double[] frequencies);

            void setupDifferentialRates(double[] differentialRates, double[] relativeRates, double normalizingConstant);
        }
    }
}
