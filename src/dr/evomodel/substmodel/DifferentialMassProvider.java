/*
 * DifferentialMassProvider.java
 *
 * Copyright © 2002-2024 the BEAST Development Team
 * http://beast.community/about
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
 *
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
        EXACT("exact") {
            @Override
            public double[] dispatch(double time,
                                     DifferentiableSubstitutionModel model,
                                     WrappedMatrix infinitesimalDifferentialMatrix) {

                return DifferentiableSubstitutionModelUtil.getExactDifferentialMassMatrix(
                        time, infinitesimalDifferentialMatrix, model.getEigenDecomposition());
            }

            @Override
            public String getReport() {
                return "Exact";
            }
        },
        FIRST_ORDER("firstOrder") {
            @Override
            public double[] dispatch(double time,
                                     DifferentiableSubstitutionModel model,
                                     WrappedMatrix infinitesimalDifferentialMatrix) {

                return DifferentiableSubstitutionModelUtil.getApproximateDifferentialMassMatrix(
                        time, infinitesimalDifferentialMatrix);
            }

            @Override
            public String getReport() {
                return "Approximate wrt parameter";
            }
        },
        AFFINE("affineCorrected") {
            @Override
            public double[] dispatch(double time,
                                     DifferentiableSubstitutionModel model,
                                     WrappedMatrix infinitesimalDifferentialMatrix) {

                return DifferentiableSubstitutionModelUtil.getAffineDifferentialMassMatrix(
                        time, infinitesimalDifferentialMatrix, model.getEigenDecomposition());
            }

            @Override
            public String getReport() {
                return "Affine-corrected wrt parameter";
            }
        };

        private final String label;

        Mode(String label) {
            this.label = label;
        }

        public String getLabel() {
            return label;
        }

        public abstract double[] dispatch(double time,
                                   DifferentiableSubstitutionModel model,
                                   WrappedMatrix infinitesimalDifferentialMatrix);

        public abstract String getReport();

        public static Mode parse(String label) {
            for (Mode mode : Mode.values()) {
                if (mode.label.equalsIgnoreCase(label)) {
                    return mode;
                }
            }
            throw new IllegalArgumentException("Unknown mode");
        }
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
