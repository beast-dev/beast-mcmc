/*
 * MasBirthDeathSerialSamplingModel.java
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

public class MasBirthDeathSerialSamplingModel extends NewBirthDeathSerialSamplingModel {

    public MasBirthDeathSerialSamplingModel(Parameter birthRate, Parameter deathRate, Parameter serialSamplingRate, Parameter treatmentProbability, Parameter samplingProbability, Parameter originTime, boolean condition, int numIntervals, double gridEnd, Type units) {
        super(birthRate, deathRate, serialSamplingRate, treatmentProbability, samplingProbability, originTime, condition, numIntervals, gridEnd, units);
    }

    @Override
    public final double processModelSegmentBreakPoint(int model, double intervalStart, double intervalEnd, int nLineages) {
//        double lnL = nLineages * (logQ(model, intervalEnd) - logQ(model, intervalStart));
        double lnL = nLineages * Math.log(Q(model, intervalEnd) / Q(model, intervalStart));
        if ( samplingProbability.getValue(model + 1) > 0.0 && samplingProbability.getValue(model + 1) < 1.0) {
            // Add in probability of un-sampled lineages
            // We don't need this at t=0 because all lineages in the tree are sampled
            // TODO: check if we're right about how many lineages are actually alive at this time. Are we inadvertently over-counting or under-counting due to samples added at this _exact_ time?
            lnL += nLineages * Math.log(1.0 - samplingProbability.getValue(model + 1));
        }
        this.savedLogQ = Double.NaN;
        return lnL;
    }

    void accumulateGradientForInterval(final double[] gradient, final int currentModelSegment, final int nLineages,
                                       final double[] partialQ_all_old, final double Q_Old,
                                       final double[] partialQ_all_young, final double Q_young) {

        for (int k = 0; k <= currentModelSegment; k++) {
            gradient[k * 5 + 0] += nLineages * (partialQ_all_old[k * 4 + 0] / Q_Old
                    - partialQ_all_young[k * 4 + 0] / Q_young);
            gradient[k * 5 + 1] += nLineages * (partialQ_all_old[k * 4 + 1] / Q_Old
                    - partialQ_all_young[k * 4 + 1] / Q_young);
            gradient[k * 5 + 2] += nLineages * (partialQ_all_old[k * 4 + 2] / Q_Old
                    - partialQ_all_young[k * 4 + 2] / Q_young);
        }
    }

    void accumulateGradientForSerialSampling(double[] gradient, int currentModelSegment, double term1,
                                             double[] intermediate) {

        for (int k = 0; k <= currentModelSegment; k++) {
                gradient[k * 5 + 0] += term1 * intermediate[k * 4 + 0];
                gradient[k * 5 + 1] += term1 * intermediate[k * 4 + 1];
                gradient[k * 5 + 2] += term1 * intermediate[k * 4 + 2];
        }
    }

     void accumulateGradientForIntensiveSampling(double[] gradient, int currentModelSegment, double term1,
                                                   double[] intermediate) {

        for (int k = 0; k < currentModelSegment; k++) {
                gradient[k * 5 + 0] += term1 * intermediate[k * 4 + 0];
                gradient[k * 5 + 1] += term1 * intermediate[k * 4 + 1];
                gradient[k * 5 + 2] += term1 * intermediate[k * 4 + 2];
        }
    }

     void dBCompute(int model, double[] dB) {

        for (int k = 0; k < model; ++k) {
            for (int p = 0; p < 4; p++) {
                dB[k * 4 + p] = -2 * (1 - rho) * lambda / A * dPModelEnd[k * 4 + p];
            }
        }

        double term1 = 1 - 2 * (1 - rho) * previousP;

        dB[model * 4 + 0] = (A * term1 - dA[0] * (term1 * lambda + mu + psi)) / (A * A);
        dB[model * 4 + 1] = (A - dA[1] * (term1 * lambda + mu + psi)) / (A * A);
        dB[model * 4 + 2] = (A - dA[2] * (term1 * lambda + mu + psi)) / (A * A);
    }

    void dPCompute(int model, double t, double intervalStart, double eAt, double[] dP, double[] dG2) {

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
        dP[model * 4 + 1] = (1 - dG2[1]) / (2 * lambda);
        dP[model * 4 + 2] = (1 - dG2[2]) / (2 * lambda);
    }


    void dQCompute(int model, double t, double[] dQ, double eAt) {

        double dwell = t - modelStartTimes[model];
        double G1 = g1(eAt);

        double term1 = 8 * eAt;
        double term2 = G1 / 2 - eAt * (1 + B);
        double term3 = eAt - 1;
        double term4 = G1 * G1 * G1;
        double term5 = -term1 * term3 / term4;

        for (int k = 0; k < model; ++k) {
            dQ[k * 4 + 0] = term5 * dB[k * 4 + 0];
            dQ[k * 4 + 1] = term5 * dB[k * 4 + 1];
            dQ[k * 4 + 2] = term5 * dB[k * 4 + 2];
        }

        double term6 = term1 / term4;
        double term7 = dwell * term2;

        dQ[model * 4 + 0] = term6 * (dA[0] * term7 - dB[model * 4 + 0] * term3);
        dQ[model * 4 + 1] = term6 * (dA[1] * term7 - dB[model * 4 + 1] * term3);
        dQ[model * 4 + 2] = term6 * (dA[2] * term7 - dB[model * 4 + 2] * term3);
    }


    final double Q(int model, double time) {
        double At = A * (time - modelStartTimes[model]);
        double eAt = Math.exp(At);
        double sqrtDenominator = g1(eAt);
        return eAt / (sqrtDenominator * sqrtDenominator);
    }

    final double logQ(int model, double time) {
        double At = A * (time - modelStartTimes[model]);
        double eAt = Math.exp(At);
        double sqrtDenominator = g1(eAt);
        return At - 2 * Math.log(sqrtDenominator); // TODO log4 (additive constant) is not needed since we always see logQ(a) - logQ(b)
    }

    @Override
    public double processInterval(int model, double tYoung, double tOld, int nLineages) {
        double logQ_young;
        double logQ_old = Q(model, tOld);
        if (!Double.isNaN(this.savedLogQ)) {
            logQ_young = this.savedLogQ;
        } else {
            logQ_young = Q(model, tYoung);
        }
        this.savedLogQ = logQ_old;
        return nLineages * Math.log(logQ_old / logQ_young);
    }

    @Override
    public double processSampling(int model, double tOld) {

        double logSampProb;

        boolean sampleIsAtEventTime = tOld == modelStartTimes[model];
        boolean samplesTakenAtEventTime = rho > 0;

        if (sampleIsAtEventTime && samplesTakenAtEventTime) {
            logSampProb = Math.log(rho);
            if (model > 0) {
                logSampProb += Math.log(r + ((1.0 - r) * previousP));
            }
        } else {
            double logPsi = Math.log(psi);
            logSampProb = logPsi + Math.log(r + (1.0 - r) * p(model,tOld));
        }

        return logSampProb;
    }
}

/*
 * Notes on inlining:
 *  https://www.baeldung.com/jvm-method-inlining#:~:text=Essentially%2C%20the%20JIT%20compiler%20tries,times%20we%20invoke%20the%20method.
 *  https://miuv.blog/2018/02/25/jit-optimizations-method-inlining/
 * static, private, final 
 */
