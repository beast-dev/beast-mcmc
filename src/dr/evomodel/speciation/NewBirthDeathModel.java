/*
 * CriticalBirthDeathSerialSamplingModel.java
 *
 * Copyright (c) 2002-2022 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.evomodel.speciation;

import dr.inference.model.Parameter;

import static java.lang.Math.exp;
import static java.lang.Math.log;


public class NewBirthDeathModel extends NewBirthDeathSerialSamplingModel {

    private int n_events;
    private double[] temp2;
    private double[] temp3;
    private double saved_q;
    private boolean partialqKnown;
    private double[] partialq;

    public NewBirthDeathModel(
            String modelName,
            Parameter birthRate,
            Parameter deathRate,
            Parameter serialSamplingRate,
            Parameter treatmentProbability,
            Parameter samplingFractionAtPresent,
            Parameter originTime,
            boolean condition,
            Type units) {

        super(modelName, birthRate, deathRate, serialSamplingRate, treatmentProbability, samplingFractionAtPresent, originTime, condition, units);
        n_events = 0;
        this.temp2 = new double[4];
        this.temp3 = new double[4];
        this.saved_q = Double.MIN_VALUE;
        this.partialqKnown = false;
        this.partialq = new double[4];
    }
    public double q(double t){
        double C1 = getC1();
        double C2 = getC2();
        // TODO why the factor of 4 and inversion here?
        double expC1t = Math.exp(C1 * t);
        return 4/(2.0 * (1.0 - Math.pow(C2,2.0)) + (1.0/expC1t) * Math.pow((1.0 - C2),2.0) + expC1t * Math.pow(1.0 + C2,2.0));
    }
    public double[] partialqpartialAll (double[] partialQ_all, double t) {
        double Q = Q(t);
        double[] temp1 = partialQpartialAll(partialQ_all, t);
        for (int i = 0; i < temp1.length; i++) {
            temp1[i] *= -4*Math.pow(Q, -2);
        }
        return temp1;
    }

    @Override
    public void precomputeConstants() {
        this.C1 = c1(lambda(), mu(), psi());
        this.C2 = c2(lambda(), mu(), psi(), rho());
        n_events = 0;
    }

    @Override
    public double processInterval(int model, double tYoung, double tOld, int nLineages) {
        // TODO Do something different
        return super.processInterval(model, tYoung, tOld, nLineages);
    }

    @Override
    public double processOrigin(int model, double rootAge) {
        double lambda = lambda();
        double rho = rho();
        double mu = mu();
        double v = exp(-(lambda - mu) * rootAge);
        double p_n = log(lambda*rho + (lambda*(1-rho) - mu)* v) - log(1- v);
        return -2*logq(rootAge) + (n_events-1)*p_n;
    }

    @Override
    public double processCoalescence(int model, double tOld) {
        n_events += 1;
        return 0;
    }

    @Override
    public double processSampling(int model, double tOld) {
        return 0;
    }

    @Override
    public double logConditioningProbability() {
        return -log(n_events);
    }

    @Override
    public void precomputeGradientConstants() {
        n_events = 0;
        this.saved_q = Double.MIN_VALUE;
        this.partialqKnown = false;
    }

    @Override
    public void processGradientInterval(double[] gradient, int currentModelSegment, double intervalStart, double intervalEnd, int nLineages) {
        double tOld = intervalEnd;
        double tYoung = intervalStart;
        double[] partialq_all_old = partialqpartialAll(temp2, tOld);
        double[] partialq_all_young;
        double q_Old = q(tOld);
        double q_young;
        if (this.saved_q != Double.MIN_VALUE) {
            q_young = this.saved_q;
        }
        else {
            q_young = q(tYoung);
        }
        this.saved_q = q_Old;

        if (partialqKnown) {
            partialq_all_young = temp3;
            System.arraycopy(partialq, 0, partialq_all_young, 0, 4);
        } else {
            partialq_all_young = partialQpartialAll(temp3, tYoung);
            //System.arraycopy(partialQ_all_young, 0, savedPartialQ, 0, 4);
            partialqKnown = true;
        }
        System.arraycopy(partialq_all_old, 0, partialq, 0, 4);


        for (int j = 0; j < 4; ++j) {
            gradient[j] += nLineages*(partialq_all_young[j] / q_young - partialq_all_old[j] / q_Old);
        }

    }

    @Override
    public void processGradientSampling(double[] gradient, int currentModelSegment, double intervalEnd) {
        return;
    }

    @Override
    public void processGradientCoalescence(double[] gradient, int currentModelSegment, double intervalEnd) {
        n_events += 1;
        return;
    }

    @Override
    public void processGradientOrigin(double[] gradient, int currentModelSegment, double totalDuration) {
        double lambda = lambda();
        double rho = rho();
        double mu = mu();
        double v = exp(-(lambda - mu) * totalDuration);
        double v2 = (lambda*(1-rho) - mu)* v;
        double v1 = lambda*rho + v2;

        // (lambda, mu, psi, rho)
        gradient[0] += (n_events - 1)*((1/v1)*(rho + v*(1-rho) - v2*totalDuration) - 1/(1-v)*v*totalDuration);
        gradient[1] +=  (n_events - 1)*(1/v1*(-v+v2*totalDuration) + 1/(1-v)*v*totalDuration);
        gradient[3] += (n_events - 1)*(1/v1*(lambda-v*lambda));

        double[] partialq_all_root = partialqpartialAll(temp3, totalDuration);

        double q_totalDuration = q(totalDuration);
        for (int i = 0; i < 4; ++i) {
            gradient[i] -=  2* partialq_all_root[i] / q_totalDuration;
        }
    }

    @Override
    public void logConditioningProbability(double[] gradient) {
        return;
    }

}
