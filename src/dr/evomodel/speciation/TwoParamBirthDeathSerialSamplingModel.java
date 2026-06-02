/*
 * TwoParamBirthDeathSerialSamplingModel.java
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

package dr.evomodel.speciation;

import dr.inference.model.Parameter;

public class TwoParamBirthDeathSerialSamplingModel extends MasBirthDeathSerialSamplingModel {

    public TwoParamBirthDeathSerialSamplingModel(Parameter birthRate, Parameter deathRate, Parameter serialSamplingRate, Parameter treatmentProbability, Parameter samplingProbability, Parameter originTime, boolean condition, int numIntervals, double gridEnd, Type units) {
        super(birthRate, deathRate, serialSamplingRate, treatmentProbability, samplingProbability, originTime, condition, numIntervals, gridEnd, units);
    }


    final void accumulateGradientForInterval(final double[] gradient, final int currentModelSegment, final int nLineages,
                                             final double[] partialQ_all_old, final double Q_Old,
                                             final double[] partialQ_all_young, final double Q_young) {

        for (int k = 0; k <= currentModelSegment; k++) {
            gradient[k * 5 + 0] += nLineages * (partialQ_all_old[k * 4 + 0] / Q_Old
                    - partialQ_all_young[k * 4 + 0] / Q_young);
            gradient[k * 5 + 2] += nLineages * (partialQ_all_old[k * 4 + 2] / Q_Old
                    - partialQ_all_young[k * 4 + 2] / Q_young);
        }
    }

    final void accumulateGradientForSerialSampling(double[] gradient, int currentModelSegment, double term1,
                                             double[] intermediate) {

        for (int k = 0; k <= currentModelSegment; k++) {
                gradient[k * 5 + 0] += term1 * intermediate[k * 4 + 0];
                gradient[k * 5 + 2] += term1 * intermediate[k * 4 + 2];
        }

    }

    final void accumulateGradientForIntensiveSampling(double[] gradient, int currentModelSegment, double term1,
                                                   double[] intermediate) {

        for (int k = 0; k < currentModelSegment; k++) {
                gradient[k * 5 + 0] += term1 * intermediate[k * 4 + 0];
                gradient[k * 5 + 2] += term1 * intermediate[k * 4 + 2];
        }

    }

    final void dBCompute(int model, double[] dB) {

        for (int k = 0; k < model; ++k) {
            for (int p = 0; p < 4; p++) {
                dB[k * 4 + p] = -2 * (1 - rho) * lambda / A * dPModelEnd[k * 4 + p];
            }
        }

        double term1 = 1 - 2 * (1 - rho) * previousP;

        dB[model * 4 + 0] = (A * term1 - dA[0] * (term1 * lambda + mu + psi)) / (A * A);
        dB[model * 4 + 2] = (A - dA[2] * (term1 * lambda + mu + psi)) / (A * A);
    }

    final void dPCompute(int model, double t, double intervalStart, double eAt, double[] dP, double[] dG2) {

        double G1 = g1(eAt);

        double term1 = -A / lambda * ((1 - B) * (eAt - 1) + G1) / (G1 * G1);

        for (int k = 0; k  < model; k ++) {
            for (int p = 0; p < 4; p++) {
                dP[k * 4 + p] = term1 * dB[k * 4 + p];
            }
        }

        for (int p = 0; p < 3; ++p) {
            double term2 = eAt * (1 + B) * dA[p] * (t - intervalStart) + (eAt - 1) * dB[model * 4 + p];
            dG2[p] = dA[p] - 2 * (G1 * (dA[p] * (1 - B) - dB[model * 4 + p] * A) - (1 - B) * term2 * A) / (G1 * G1);
        }

        double G2 = g2(G1);

        dP[model * 4 + 0] = (-mu - psi - lambda * dG2[0] + G2) / (2 * lambda * lambda);
        dP[model * 4 + 2] = (1 - dG2[2]) / (2 * lambda);
    }


    final void dQCompute(int model, double t, double[] dQ, double eAt) {

        double dwell = t - modelStartTimes[model];
        double G1 = g1(eAt);

        double term1 = 8 * eAt;
        double term2 = G1 / 2 - eAt * (1 + B);
        double term3 = eAt - 1;
        double term4 = G1 * G1 * G1;
        double term5 = -term1 * term3 / term4;

        for (int k = 0; k < model; ++k) {
            dQ[k * 4 + 0] = term5 * dB[k * 4 + 0];
            dQ[k * 4 + 2] = term5 * dB[k * 4 + 2];
        }

        double term6 = term1 / term4;
        double term7 = dwell * term2;

        dQ[model * 4 + 0] = term6 * (dA[0] * term7 - dB[model * 4 + 0] * term3);
        dQ[model * 4 + 2] = term6 * (dA[2] * term7 - dB[model * 4 + 2] * term3);
    }

}

/*
 * Notes on inlining:
 *  https://www.baeldung.com/jvm-method-inlining#:~:text=Essentially%2C%20the%20JIT%20compiler%20tries,times%20we%20invoke%20the%20method.
 *  https://miuv.blog/2018/02/25/jit-optimizations-method-inlining/
 * static, private, final 
 */
