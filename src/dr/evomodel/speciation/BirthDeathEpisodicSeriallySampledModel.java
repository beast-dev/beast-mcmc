/*
 * BirthDeathSerialSkylineModel.java
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

package dr.evomodel.speciation;

import dr.evolution.io.Importer;
import dr.evolution.io.NewickImporter;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.util.Taxon;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;
import dr.util.Citable;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Set;

/**
 * The episodic serially-sampled birth-death model of Gavryushkina et al. (2014)
 * If sampled ancestors are ever added, this will need overhauling!
 *
 * @author Andy Magee
 */
public class BirthDeathEpisodicSeriallySampledModel extends SpeciationModel implements SpeciationModelGradientProvider{


    // extant sampling proportion
    Parameter samplingFractionAtPresent;

    // birth rate
    Parameter birthRate;

    // death rate
    Parameter deathRate;

    // serial sampling rate
    Parameter serialSamplingRate;

    // sampling probability at event times
    Parameter samplingProbability;

    // "treatmentProbability" parameter aka r aka Pr(death | lineage is sampled)
    Parameter treatmentProbability;

    // the originTime of the infection, origin > tree.getRoot();
    Parameter originTime;

    // TODO if we want to supplant other birth-death models, need an ENUM, and choice of options
    // Minimally, need survival of 1 lineage (passable default for SSBDP) and nTaxa (which is current option for non-serially-sampled BDP)
    private boolean conditionOnSurvival;

    private boolean birthRateChanges = false;
    private boolean deathRateChanges = false;
    private boolean serialSamplingRateChanges = false;
    private boolean treatmentChanges = false;
    private boolean intensiveSamplingOnlyAtPresent = false;
    private boolean noIntensiveSampling = true;


    boolean computedBCurrent;
    private double[][] partialBCurrentPartialAll;

    private double[][] partialPPreviousPartialAll;

    private double[][] partialPCurrentPartialAll;


    // Tolerance for declaring that a node time is equal to an event time
    double absTol = 1e-8;

    int numIntervals = 1;
    double gridEnd;
    // there are numIntervals intervalTimes, implicitly intervalTimes[-1] == 0
    double[] intervalTimes;

    // Pre-computed values which are functions of parameters
    // TODO precompute q_{i-1}(ti)
    protected double[] piMinus1;
    protected double[] Ai;
    protected double[] Bi;

    public BirthDeathEpisodicSeriallySampledModel(
            Parameter birthRate,
            Parameter deathRate,
            Parameter serialSamplingRate,
            Parameter r,
            Parameter samplingProbability,
            Parameter originTime,
            int numIntervals,
            double gridEnd,
            Type units) {

        this("birthDeathEpisodicSeriallySampledModel", birthRate, deathRate, serialSamplingRate, r, samplingProbability,
                originTime, numIntervals, gridEnd, units);
    }

    public SpeciationModelGradientProvider getProvider() { // This is less INTRUSIVE to the exisiting file
        return this;
    }

    public BirthDeathEpisodicSeriallySampledModel(
            String modelName,
            Parameter birthRate,
            Parameter deathRate,
            Parameter serialSamplingRate,
            Parameter treatmentProbability,
            Parameter samplingProbability,
            Parameter originTime,
            int numIntervals,
            double gridEnd,
            Type units) {

        super(modelName, units);

        this.numIntervals = numIntervals;
        this.gridEnd = gridEnd;
        setupTimeline();

        if (birthRate.getSize() != 1 && birthRate.getSize() != numIntervals) {
            throw new RuntimeException("Length of birthRate parameter should be one or equal to the size of time parameter (size = " + numIntervals + ")");
        }

        if (deathRate.getSize() != 1 && deathRate.getSize() != numIntervals) {
            throw new RuntimeException("Length of deathRate parameter should be one or equal to the size of time parameter (size = " + numIntervals + ")");
        }

        if (serialSamplingRate.getSize() != 1 && serialSamplingRate.getSize() != numIntervals) {
            throw new RuntimeException("Length of serialSamplingRate parameter should be one or equal to the size of time parameter (size = " + numIntervals + ")");
        }

        if (treatmentProbability.getSize() != 1 && treatmentProbability.getSize() != numIntervals) {
            throw new RuntimeException("Length of r parameter should be one or equal to the size of time parameter (size = " + numIntervals + ")");
        }

        if (samplingProbability.getSize() != 1 && samplingProbability.getSize() != (numIntervals)) {
            throw new RuntimeException("Length of samplingProbability parameter should be one or equal to the size of time parameter (size = " + numIntervals + ")");
        }

        if (birthRate.getSize() > 1) {
            birthRateChanges = true;
        }

        if (deathRate.getSize() > 1) {
            deathRateChanges = true;
        }

        if (serialSamplingRate.getSize() > 1) {
            serialSamplingRateChanges = true;
        }

        if (treatmentProbability.getSize() > 1) {
            treatmentChanges = true;
        }

        if (samplingProbability.getSize() > 1) {
            intensiveSamplingOnlyAtPresent = false;
        }

        for (int i=0; i < samplingProbability.getSize(); i++) {
            if ( samplingProbability.getValue(i) > Double.MIN_VALUE ) {
                noIntensiveSampling = false;
                break;
            }
        }

        this.birthRate = birthRate;
        addVariable(birthRate);
        birthRate.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, birthRate.getSize()));

        this.deathRate = deathRate;
        addVariable(deathRate);
        deathRate.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, deathRate.getSize()));

        this.serialSamplingRate = serialSamplingRate;
        addVariable(serialSamplingRate);
        deathRate.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, serialSamplingRate.getSize()));

        this.treatmentProbability = treatmentProbability;
        addVariable(treatmentProbability);
        samplingProbability.addBounds(new Parameter.DefaultBounds(1.0, 0.0, treatmentProbability.getSize()));

        this.samplingProbability = samplingProbability;
        addVariable(samplingProbability);
        samplingProbability.addBounds(new Parameter.DefaultBounds(1.0, 0.0, samplingProbability.getSize()));

        this.originTime = originTime;
        addVariable(originTime);
        originTime.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, originTime.getSize()));

    }

    // TODO should probably be replaced and brought in line with smoothSkygrid
    private void setupTimeline() {
        if (intervalTimes == null) {
            intervalTimes = new double[numIntervals];
        } else {
            Arrays.fill(intervalTimes, 0.0);
        }

        for (int idx = 0; idx <= numIntervals - 1 ; idx++) {
            intervalTimes[idx] = (idx+1) * (gridEnd / numIntervals);
        }
    }

    public double Ai(double lambda, double mu, double psi) {

        return Math.sqrt(Math.pow((lambda - mu - psi),2.0) + 4.0 * lambda * psi);
    }

    public double Bi(double lambda, double mu, double psi, double rho, double A, double pPrevious) {

        return ((1.0 - 2.0 * (1.0 - rho)*pPrevious) * lambda + mu + psi)/A;
    }

    public double p(int index, double t) {
        double lambda = birthRate(index);
        double mu = deathRate(index);
        double psi = serialSamplingRate(index);
        double A = Ai[index];
        double B = Bi[index];

        double ti = index == 0 ? 0 : intervalTimes[index-1];
        double eA = Math.exp(A * (t - ti));
        double oneMinus = eA * (1.0 + B) - (1.0 - B);
        double onePlus = eA * (1.0 + B) + (1.0 - B);
        //System.out.println((lambda + mu + psi- A * oneMinus/onePlus)/(2 * lambda));
        return (lambda + mu + psi- A * oneMinus/onePlus)/(2 * lambda);

    }

    // interval i is (t_{i-1},t_i], except for the first interval which also includes 0
    public double p(int index, double t, double lambda, double mu, double rho, double A, double B) {
        double ti = index == 0 ? 0 : intervalTimes[index-1];
        double eA = Math.exp(A * (t - ti));
        double oneMinus = eA * (1.0 + B) - (1.0 - B);
        double onePlus = eA * (1.0 + B) + (1.0 - B);
        //System.out.println((lambda + mu + rho - A * oneMinus/onePlus)/(2 * lambda));
        return (lambda + mu + rho - A * oneMinus/onePlus)/(2 * lambda);

    }

    public double q(int index, double t) {
        double ti = index == 0 ? 0 : intervalTimes[index-1];

        double eA = Math.exp(Ai[index] * (t - ti));
        double sqrtDenom = eA * (1.0 + Bi[index]) + (1.0 - Bi[index]);

        return (4.0 * eA)/(Math.pow(sqrtDenom,2.0));
    }

    public double logq(int index, double t) {
        double ti = index == 0 ? 0 : intervalTimes[index-1];

        double At = Ai[index] * (t - ti);
        double eA = Math.exp(At);
        double sqrtDenom = eA * (1.0 + Bi[index]) + (1.0 - Bi[index]);
        //System.out.println(At + Math.log(4.0) - 2.0 * Math.log(sqrtDenom));
        return At + Math.log(4.0) - 2.0 * Math.log(sqrtDenom);
    }

    public double birthRate(int i) {
        return birthRate.getValue(birthRateChanges ? i : 0);
    }

    public double deathRate(int i) {
        return deathRate.getValue(deathRateChanges ? i : 0);
    }

    public double serialSamplingRate(int i) {
        return serialSamplingRate.getValue(serialSamplingRateChanges ? i : 0);
    }

    public double r(int i) {
        return treatmentProbability.getValue(treatmentChanges ? i : 0);
    }

    public double samplingProbability(int i) {
        if ( intensiveSamplingOnlyAtPresent ) {
            if ( i == 0) {
                return samplingProbability.getValue(0);
            } else {
                return 0.0;
            }
        } else {
            return samplingProbability.getValue(i);
        }
    }


    private double lambda(int i) {
        return birthRate.getValue(birthRateChanges ? i : 0);
    }

    private double mu(int i) {
        return deathRate.getValue(deathRateChanges ? i : 0);
    }

    private double psi(int i) {
        return serialSamplingRate.getValue(serialSamplingRateChanges ? i : 0);
    }

    private double rho(int i) {
        if ( intensiveSamplingOnlyAtPresent ) {
            if ( i == 0) {
                return samplingProbability.getValue(0);
            } else {
                return 0.0;
            }
        } else {
            return samplingProbability.getValue(i);
        }
    }

    /*    calculate and store A_i, B_i and p_i(t_i)        */
    public void updateModelValues(int model) {

        double t_origin = originTime.getValue(0);

        Ai = new double[numIntervals];
        Bi = new double[numIntervals];
        piMinus1 = new double[numIntervals];

        for (int i = 0; i < numIntervals; i++) {
            Ai[i] = Ai(birthRate(i), deathRate(i), serialSamplingRate(i));
            //System.out.println("Ai[" + i + "]=" + Ai[i]);
        }

        piMinus1[0] = 1.0;
        Bi[0] = Bi(birthRate(0), deathRate(0), serialSamplingRate(0), samplingProbability(0), Ai[0], piMinus1[0]);
        //System.out.println("Bi[0]=" + Bi[0]);
        for (int i = 1; i < numIntervals; i++) {
            piMinus1[i] = p(i-1, intervalTimes[i-1]);
            Bi[i] = Bi(birthRate(i), deathRate(i), serialSamplingRate(i), samplingProbability(i), Ai[i], piMinus1[i]);
        }

        if (partialBCurrentPartialAll == null) {
            partialPCurrentPartialAll = new double[numIntervals][4];
            partialBCurrentPartialAll = new double[numIntervals][4];
        }

        computedBCurrent = false;
        // computedPPrevious = false;
        //partialBPreviousPartialAll = partialBCurrentPartialAll;
        partialPPreviousPartialAll = partialPCurrentPartialAll;
        // Use deep copy here, better to transform it into 1d array;
//            for (int k=0; k<model; k++) {
//                for (int n=0; n<4; n++) {
//                    partialPPreviousPartialAll[k][n] = partialPCurrentPartialAll[k][n];
//                }
//            }
        for (int k = 0; k <= model; k++) {
            partialPCurrentPartialAll[k] = partialPpartialAll(model, k, intervalTimes[model]);
        }


    }

    /**
     * Likelihood calculation
     *
     * @param tree the tree to calculate likelihood of
     * @return log-likelihood of density
     */
    public final double calculateTreeLogLikelihood(Tree tree) {
        throw new RuntimeException("Not yet implemented!");
    }

    public double calculateTreeLogLikelihood(Tree tree, Set<Taxon> exclude) {
        if (exclude.size() == 0) return calculateTreeLogLikelihood(tree);
        throw new RuntimeException("Not implemented!");
    }


    @Override
    public double[] getBreakPoints() {
        return intervalTimes;
    }

    @Override
    public double processModelSegmentBreakPoint(int model, double intervalStart, double segmentIntervalEnd, int nLineages) {
        return nLineages * (logq(model, segmentIntervalEnd) - logq(model, intervalStart));
    }

    @Override
    public double processInterval(int model, double tYoung, double tOld, int nLineages) {
        return nLineages * (logq(model, tOld) - logq(model, tYoung));
    }

    @Override
    public double processOrigin(int model, double rootAge) {
        return (logq(model,originTime.getValue(0)) - logq(model, rootAge));
    }

    @Override
    public double processCoalescence(int model, double tOld) {
        return Math.log(birthRate.getValue(model)); // TODO Notice the natural parameterization is `log lambda`
    }


    @Override
    public double processSampling(int model, double tOld) {
        // TODO some of these might work better as stored variables to avoid recomputation
        double logSampProb = 0.0;

        double atEventTimeTolerance = Double.MIN_VALUE;

        boolean sampleIsAtPresent = tOld < atEventTimeTolerance;
        boolean samplesTakenAtPresent = samplingProbability(0) >= Double.MIN_VALUE;

        boolean sampleIsAtEventTime = Math.abs(tOld - intervalTimes[model]) < atEventTimeTolerance;
        boolean samplesTakenAtEventTime = samplingProbability(model) >= Double.MIN_VALUE;

        if (sampleIsAtPresent && samplesTakenAtPresent) {
            logSampProb = Math.log(samplingProbability(0));
        } else if (sampleIsAtEventTime && samplesTakenAtEventTime) {
            logSampProb = Math.log(samplingProbability(model+1));
        } else {
            double logPsi = Math.log(serialSamplingRate(model)); // TODO Notice the natural parameterization is `log psi`
            double r = treatmentProbability.getValue(model);
            logSampProb = logPsi + Math.log(r + (1.0 - r) * p(model,tOld));
        }

        return logSampProb;
    }

    private double ti(int i) {
        return i == 0 ? 0 : intervalTimes[i-1];
    }

    @Override
    public double logConditioningProbability() {
       return 0;
    }

    private double partialApartialLambda(int i) {
        return (lambda(i) - mu(i) + psi(i)) / Ai[i];

    }

    private double partialApartialMu(int i) {
        return (-lambda(i) + mu(i) + psi(i)) / Ai[i];
    }

    private double partialApartialPsi(int i) {
        return (lambda(i) + mu(i) + psi(i)) / Ai[i];
    }

    private double[] partialApartialAll(int i) {
        double[] partialA = new double[4];
        partialA[0] = partialApartialLambda(i);
        partialA[1] = partialApartialMu(i);
        partialA[2] = partialApartialPsi(i);
        partialA[3] = 0;
        return partialA;
    }


    private double[] partialBpartialAll(int i, int k) {
        if (computedBCurrent) {
            return partialBCurrentPartialAll[k];
        }
        double[] partialB = new double[4];
        if (k == i) {
            double[] partialA = partialApartialAll(i);
            double temp = 1 - 2 * (1 - rho(i)) * piMinus1[i];
            partialB[0] = (Ai[i] * temp - partialA[0] * (temp * lambda(i) + mu(i) + psi(i))) / (Ai[i] * Ai[i]);
            partialB[1] = (Ai[i] - partialA[1] * (temp * lambda(i) + mu(i) + psi(i))) / (Ai[i] * Ai[i]);
            partialB[2] = (Ai[i] - partialA[2] * (temp * lambda(i) + mu(i) + psi(i))) / (Ai[i] * Ai[i]);
            partialB[3] = 2 * lambda(i) * piMinus1[i] / Ai[i];
        } else if (k < i) {
            for(int n = 0; n < 4; n++) {
                partialB[n] = -2 * (1 - rho(i)) * lambda(i) / Ai[i] * partialPPreviousPartialAll[k][n]; // partialPpartialAll(i-1, k, ti, true)[n];
            }
        }
        // else: k > i, all zero
        partialBCurrentPartialAll[k] = partialB;
        if (k == i) {
            computedBCurrent = true;
        }
        return partialB;
    }


    private double[] partialPpartialAll(int i, int k, double t) {
//        if (computedPPrevious && save) {
//            return partialPPreviousPartialAll[k];
//        }
        double[] partialP = new double[4];
        if (k == i) {
            double ti = i == 0 ? 0 : intervalTimes[i-1];
            double g1 = Math.exp((Ai[i]) * (t - ti)) * (1 + Bi[i]) + (1 - Bi[i]);
            double g2 =  Ai[i] * (1 - 2*(1 - Bi[i]) / g1);
            double[] partialA = partialApartialAll(i);
            double[] partialB = partialBpartialAll(i, k); //partialBPreviousPartialAll[k];
            double[] temp3 = new double[3];
            for(int n = 0; n < 3; n++) {
                temp3[n] = Math.exp(Ai[i] * (t - ti)) * (1 + Bi[i]) * partialA[n] * (t - ti) + (Math.exp(Ai[i] * (t - ti)) - 1) * partialB[n];
            }
            double partialTemp1;
            // lambda
            partialTemp1 = partialA[0] - 2 * (g1 * (partialA[0] * (1 - Bi[i]) - partialB[0] * Ai[i]) -  (1 - Bi[i]) * temp3[0] * Ai[i]) / (g1 * g1);
            partialP[0] = (- mu(i) - psi(i) - lambda(i) * partialTemp1 + g2)/ (2 * lambda(i) * lambda(i));
            // mu
            partialTemp1 = partialA[1] - 2 * (g1 * (partialA[1] * (1 - Bi[i]) - partialB[1] * Ai[i]) -  (1 - Bi[i]) * temp3[1] * Ai[i]) / (g1 * g1);
            partialP[1] = (1 - partialTemp1)/(2 * lambda(i));
            // psi
            partialTemp1 = partialA[2] - 2 * (g1 * (partialA[2] * (1 - Bi[i]) - partialB[2] * Ai[i]) -  (1 - Bi[i]) * temp3[2] * Ai[i]) / (g1 * g1);
            partialP[2] = (1 - partialTemp1)/(2 * lambda(i));
            // rho
            partialP[3] = -Ai[i] / lambda(i) * ((1 - Bi[i]) * (Math.exp(Ai[i] * (t - ti)) - 1) + g1 ) * partialB[3] / Math.pow(g1,2);
        } else if (k < i) {
            double ti = i == 0 ? 0 : intervalTimes[i-1];
            double[] partialB = partialBpartialAll(i, k); // partialBPreviousPartialAll[k]; // if computed already, directly return partialBCurrentPartialAll[k]
            double g1 = Math.exp((Ai[i]) * (t - ti)) * (1 + Bi[i]) + (1 - Bi[i]);
            for (int n = 0; n < 4; n++) {
                partialP[n] = -Ai[i] / lambda(i) * ((1 - Bi[i]) * (Math.exp(Ai[i] * (t - ti)) - 1) + g1 ) * partialB[n] / Math.pow(g1,2);
            }
        }
//        if (save) {
//            partialPPreviousPartialAll[k] = partialP;
//            if (k == i) {
//                computedPPrevious = true;
//            }
//        }
        return partialP;
    }


    private double[] partialqpartialAll(int i, int k, double t) {
        // TODO: share temp2 with partialPpartialAll
        double[] partialq = new double[4];
        if (i < k) { return partialq; }
        double[] partialA = partialApartialAll(i);
        double tempA1;
        double tempA2;
        double[] partialB = partialBpartialAll(i, k);
        double ti = i == 0 ? 0 : intervalTimes[i-1];
        double g1 = Math.exp((Ai[i]) * (t - ti)) * (1 + Bi[i]) + (1 - Bi[i]);
        double temp_exp = Math.exp(Ai[i] * (t - ti));
        for (int n = 0; n < 3; n++) {
            if (i == k) {
                tempA1 = (t - ti) * partialA[n] * (g1 / 2 - temp_exp * (1 + Bi[i]));
                //tempA2 = (t - ti) * partialA[n] * temp_exp * (1 + Bi[i]);
            } else {
                tempA1 = 0;
                //tempA2 = 0;
            }
            partialq[n] = 8 * temp_exp * (tempA1 - partialB[n] * (temp_exp - 1)) / (Math.pow(g1, 3));//(tempA1 - 8 * temp_exp * temp2 * (tempA2 + partialB[n] * (temp_exp - 1))) / (Math.pow(temp2, 4));
        }
        partialq[3] = -8 * temp_exp * partialB[3] * (temp_exp - 1) / (Math.pow(g1, 3));
        return partialq;
    }

    @Override
    public Parameter getSamplingProbabilityParameter() {
        return samplingProbability;
    }

    @Override
    public Parameter getDeathRateParameter() {
        return deathRate;
    }


    @Override
    public Parameter getBirthRateParameter() {
        return birthRate;
    }



    @Override
    public Parameter getSamplingRateParameter() {
        return serialSamplingRate;
    }



    @Override
    public Parameter getTreatmentProbabilityParameter() {
        return treatmentProbability;
    }


    @Override
    public void precomputeGradientConstants() {
        //partialBPreviousPartialAll = new double[numIntervals][4];
        partialPPreviousPartialAll = new double[numIntervals][4];
        partialPCurrentPartialAll = new double[numIntervals][4];
        partialBCurrentPartialAll = new double[numIntervals][4];
        updateModelValues(0);
    }

    @Override
    public void processGradientModelSegmentBreakPoint(double[] gradient, int currentModelSegment, double intervalStart, double segmentIntervalEnd, int nLineages) {
        // only partialq_ii has value, so we don't need to loop k
        double[] partialqStart;
        double[] partialqEnd;
        double qStart;
        double qEnd;
        for (int k = 0; k <= currentModelSegment; k++) {
            partialqStart = partialqpartialAll(currentModelSegment, k, intervalStart);
            partialqEnd = partialqpartialAll(currentModelSegment, k, segmentIntervalEnd);
            qStart = q(currentModelSegment, intervalStart);
            qEnd = q(currentModelSegment, segmentIntervalEnd);

            for (int n = 0; n < 4; n++) {
                // assume lambda_1, lambda_2, ..., mu_1, mu_2, ...
                gradient[k * 5 + n] += nLineages * (partialqEnd[n] / qEnd - partialqStart[n] / qStart);
                //gradient[n * numIntervals + k] += nLineages * (partialqEnd[n] / qEnd - partialqStart[n] / qStart);
            }
        }
        // return nLineages * (logq(model, segmentIntervalEnd) - logq(model, intervalStart));
    }

    @Override
    public void processGradientInterval(double[] gradient, int currentModelSegment, double intervalStart, double intervalEnd, int nLineages) {
        // TODO: cache these calculations since they are used by processGradientModelSegmentBreakPoint
        double[] partialqStart;
        double[] partialqEnd;
        double qStart;
        double qEnd;
        for (int k = 0; k <= currentModelSegment; k++) {
            partialqStart = partialqpartialAll(currentModelSegment, k, intervalStart);
            partialqEnd = partialqpartialAll(currentModelSegment, k, intervalEnd);
            qStart = q(currentModelSegment, intervalStart);
            qEnd = q(currentModelSegment, intervalEnd);

            for (int n = 0; n < 4; n++) {
                // assume lambda_1, lambda_2, ..., mu_1, mu_2, ...
                gradient[k * 5 + n] += nLineages * (partialqEnd[n] / qEnd - partialqStart[n] / qStart);
                //gradient[n * numIntervals + k] += nLineages * (partialqEnd[n] / qEnd - partialqStart[n] / qStart);
            }
        }
        // return nLineages * (logq(model, tOld) - logq(model, tYoung));
    }

    @Override
    public void processGradientOrigin(double[] gradient, int currentModelSegment, double totalDuration) {
        double origin = originTime.getValue(0);
        double[] partialqOrigin;
        double[] partialqRoot;
        double qOrigin;
        double qRoot;

        for (int k = 0; k <= currentModelSegment; k++) {
            partialqOrigin = partialqpartialAll(currentModelSegment, k, origin);
            partialqRoot = partialqpartialAll(currentModelSegment, k, totalDuration);
            qOrigin = q(currentModelSegment, origin);
            qRoot = q(currentModelSegment, totalDuration);

            for (int n = 0; n < 4; n++) {
                // assume lambda_1, lambda_2, ..., mu_1, mu_2, ...
                gradient[k * 5 + n] += (partialqOrigin[n] / qOrigin - partialqRoot[n] / qRoot);
                //gradient[n * numIntervals + k] += (partialqOrigin[n] / qOrigin - partialqRoot[n] / qRoot);
            }
        }
        // return (logq(model,originTime.getValue(0)) - logq(model, rootAge));
    }

    @Override
    public void processGradientCoalescence(double[] gradient, int currentModelSegment, double intervalEnd) {
        gradient[currentModelSegment*5] += 1 / lambda(currentModelSegment);
        //gradient[currentModelSegment] += 1 / lambda(currentModelSegment);
        // return Math.log(birthRate.getValue(model)); // TODO Notice the natural parameterization is `log lambda`
    }


    @Override
    public void processGradientSampling(double[] gradient, int currentModelSegment, double intervalEnd) {
        // TODO some of these might work better as stored variables to avoid recomputation
        double logSampProb = 0.0;

        double atEventTimeTolerance = Double.MIN_VALUE;

        boolean sampleIsAtPresent = intervalEnd < atEventTimeTolerance;
        boolean samplesTakenAtPresent = samplingProbability(0) >= Double.MIN_VALUE;

        boolean sampleIsAtEventTime = Math.abs(intervalEnd - intervalTimes[currentModelSegment]) < atEventTimeTolerance;
        boolean samplesTakenAtEventTime = samplingProbability(currentModelSegment) >= Double.MIN_VALUE;

        if (sampleIsAtPresent && samplesTakenAtPresent) {
            // logSampProb = Math.log(samplingProbability(0));
            gradient[3] += 1 / rho(0);
            //gradient[3*numIntervals] += 1 / rho(0);
        } else if (sampleIsAtEventTime && samplesTakenAtEventTime) {
            // logSampProb = Math.log(samplingProbability(currentModelSegment+1));
            gradient[3 + 5*(currentModelSegment+1)] += 1 / rho(currentModelSegment+1);
            //gradient[3*numIntervals + currentModelSegment+1] += 1 / rho(currentModelSegment+1);
        } else {
            // double logPsi = Math.log(serialSamplingRate(currentModelSegment)); // TODO Notice the natural parameterization is `log psi`
             gradient[2 + 5*currentModelSegment] += 1 / psi(currentModelSegment);
            //gradient[2*numIntervals + currentModelSegment] += 1 / psi(currentModelSegment);
            double r = treatmentProbability.getValue(currentModelSegment);
            // logSampProb = logPsi + Math.log(r + (1.0 - r) * p(currentModelSegment,intervalEnd));
            // partialP
            double p_it = p(currentModelSegment, intervalEnd);
            gradient[4 + 5*currentModelSegment] += (1 - p_it) / ((1-r)*p_it + r);
            //gradient[4*numIntervals + currentModelSegment] += (1 - p_it) / ((1-r)*p_it + r);
            double[] partialP;
            for(int k = 0; k <= currentModelSegment; k++) {
                partialP = partialPpartialAll(currentModelSegment, k, intervalEnd);
                for(int n = 0; n < 4; n++) {
                    gradient[n + 5*k] += (1 - r) / ((1 - r) * p_it + r) * partialP[n];
                    //gradient[n*numIntervals + k] += (1 - r) / ((1 - r) * p_it + r) * partialP[n];
                }
            }
        }
    }

    @Override
    public void logConditioningProbability(double[] gradient) {
        return;
    }
    @Override
    public int getGradientLength() { return 5*numIntervals; }



}