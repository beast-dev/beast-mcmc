/*
 * NewBirthDeathSerialSamplingModel.java
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

import dr.evolution.tree.Tree;
import dr.evolution.util.Taxon;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;
import dr.util.Author;
import dr.util.Citable;
import dr.util.Citation;

import java.util.*;

/**
 * A phylogenetic birth-death-sampling model which includes serial sampling, sampling at present, and the possibility of treatmentProbability.
 */
public class NewBirthDeathSerialSamplingModel extends SpeciationModel implements SpeciationModelGradientProvider, Citable {
    // TODO should we pre-emptively put the combinatorial constant in here? It only really matters when inferring a tree where there may/may not be sampled ancestors
    // TODO don't check if rho >= minVal, just 0 (likelihood won't blow up unless it's _exactly_ 0
    // TODO the "is time t event time ti" can/should probably get moved to the point at which we pass the tree to the birth-death model
    // extant sampling proportion
    Parameter samplingProbability;

    // birth rate
    Parameter birthRate;

    // death rate
    Parameter deathRate;

    // serial sampling rate
    Parameter serialSamplingRate;

    // "treatmentProbability" parameter aka r aka Pr(death | lineage is sampled)
    Parameter treatmentProbability;

    // the originTime of the infection, origin > tree.getRoot();
    Parameter originTime;

    // Tolerance for declaring that a node time is equal to an event time
//    double absTol = 1e-8;

    // there are numIntervals intervalStarts, implicitly intervalStarts[-1] == 0
    int numIntervals;
    double gridEnd;
    double[] modelStartTimes = null;

    // TODO if we want to supplant other birth-death models, need an ENUM, and choice of options
    // Minimally, need survival of 1 lineage (passable default for SSBDP) and nTaxa (which is current option for non-serially-sampled BDP)
    private final boolean conditionOnSurvival;

    double modelStartTime;
    double A;
    double B;
    double previousP;

    double lambda;
    double mu;
    double psi;
    double r;
    double rho;
    //double rho0; // TODO remove

    private double savedQ;
    private double[] partialQ;
    private boolean partialQKnown;

    private final double[] dQStart;
    private final double[] dQEnd;

    private double[] temp33;
    private final double[] temp44;

    private double eAt_Old;
    private double eAt_End;

    private final double[] dA;
    private final double[] dB;
    private final double[] dG2;

    boolean computedBCurrent;

    private final double[] dPIntervalEnd;
    private final double[] dPModelEnd;

    public NewBirthDeathSerialSamplingModel(
            Parameter birthRate,
            Parameter deathRate,
            Parameter serialSamplingRate,
            Parameter treatmentProbability,
            Parameter samplingProbability,
            Parameter originTime,
            boolean condition,
            int numIntervals,
            double gridEnd,
            Type units) {

        this("NewBirthDeathSerialSamplingModel", birthRate, deathRate, serialSamplingRate,
                treatmentProbability, samplingProbability, originTime, condition, numIntervals, gridEnd, units);
    }

    public SpeciationModelGradientProvider getProvider() { // This is less INTRUSIVE to the exisiting file
        return this;
    }

    public NewBirthDeathSerialSamplingModel(
            String modelName,
            Parameter birthRate,
            Parameter deathRate,
            Parameter serialSamplingRate,
            Parameter treatmentProbability,
            Parameter samplingProbability,
            Parameter originTime,
            boolean condition,
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

        this.conditionOnSurvival = condition;

        this.originTime = originTime;
        if (originTime != null) {
            addVariable(originTime);
            originTime.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, 1));
        }

        this.savedQ = Double.MIN_VALUE;
        this.partialQ = new double[4 * numIntervals];
        this.partialQKnown = false;

        this.dA = new double[4];
        this.dG2 = new double[3];

        this.dB = new double[numIntervals * 4];

        this.dQStart = new double[numIntervals * 4];
        this.dQEnd = new double[numIntervals * 4];

        this.temp33 = new double[numIntervals * 4];
        this.temp44 = new double[numIntervals * 4];

        this.dPIntervalEnd = new double[numIntervals * 4];
        this.dPModelEnd = new double[numIntervals * 4];

        this.gridEnd = gridEnd;
        this.numIntervals = numIntervals;
        setupTimeline();
    }

    // TODO should probably be replaced and brought in line with smoothSkygrid
    private void setupTimeline() {
        if (modelStartTimes == null) {
            modelStartTimes = new double[numIntervals];
        } else {
            Arrays.fill(modelStartTimes, 0.0);
        }

        modelStartTimes[0] = 0;
        for (int idx = 1; idx <= numIntervals - 1 ; idx++) {
            modelStartTimes[idx] = idx * (gridEnd / numIntervals);
        }
    }

    @Override
    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        // Do nothing
    }

    private double p(int model, double t) {
        double eAt = Math.exp(A * (t - modelStartTimes[model]));
        return p(eAt);
    }

    private double p(double eAt) {
        double eAt1B = eAt * (1.0 + B);
        return (lambda + mu + psi - A * ((eAt1B - (1.0 - B)) / (eAt1B + (1.0 - B)))) / (2.0 * lambda);
    }

    private double logQ(int model, double time) {
        double At = A * (time - modelStartTimes[model]);
        double eAt = Math.exp(At);
        double sqrtDenominator = g1(eAt);
        return At + log4 - 2 * Math.log(sqrtDenominator); // TODO log4 (additive constant) is not needed since we always see logQ(a) - logQ(b)
    }

    private static final double log4 = Math.log(4.0);
    
    // Named as per Gavryushkina et al 2014
    private static double computeA(double lambda, double mu, double psi) {
        return Math.abs(Math.sqrt(Math.pow(lambda - mu - psi, 2.0) + 4.0 * lambda * psi));  // TODO Why do we have Math.abs(always positive #)?
    }

    // Named as per Gavryushkina et al 2014
    private static double computeB(double lambda, double mu, double psi, double rho, double A, double pPrevious) {
        return ((1.0 - 2.0 * (1 - rho) * pPrevious) * lambda + mu + psi)/A;
    }

    @Override
    public double logConditioningProbability() {
        double logP = 0.0;
        if ( conditionOnSurvival ) {
            double origin = originTime.getParameterValue(0);
            double[] modelBreakPoints = getBreakPoints();
            int idx = modelBreakPoints.length - 2;
            double intervalStart = modelBreakPoints[idx];
            // Origin is probably near the last index, so we do a linear search forwards in time from there
            while ( origin < intervalStart) {
                --idx;
                intervalStart = modelBreakPoints[idx];
            }
            logP -= Math.log(1.0 - p(idx, origin));
        }
        return logP;
    }

    public final double calculateTreeLogLikelihood(Tree tree) {  // TODO Remove
        throw new RuntimeException("To be removed");
    }

    @Override
    public double[] getBreakPoints() {
        // TODO Fix, this is all messed up; we should hold one set of values [0, ..., \infty]
        // TODO when fixing this, also fix EfficientSpeciationLikelihood.fixTimes()
        double[] breakPoints = new double[numIntervals];
        System.arraycopy(modelStartTimes, 1, breakPoints, 0, numIntervals - 1);

//        breakPoints[numIntervals-1] = originTime.getValue(0);
        breakPoints[numIntervals - 1] = Double.POSITIVE_INFINITY;
        return breakPoints;
    }

    @Override
    public double processModelSegmentBreakPoint(int model, double intervalStart, double intervalEnd, int nLineages) {
        double lnL = nLineages * (logQ(model, intervalEnd) - logQ(model, intervalStart));
        if ( samplingProbability.getValue(model + 1) > 0.0 && samplingProbability.getValue(model + 1) < 1.0) {
            // Add in probability of un-sampled lineages
            // We don't need this at t=0 because all lineages in the tree are sampled
            // TODO: check if we're right about how many lineages are actually alive at this time. Are we inadvertently over-counting or under-counting due to samples added at this _exact_ time?
            lnL += nLineages * Math.log(1.0 - samplingProbability.getValue(model + 1));
        }
        return lnL;
    }

    private void updateParameterValues(int model) {
        lambda = birthRate.getParameterValue(model);
        mu = deathRate.getParameterValue(model);
        psi = serialSamplingRate.getParameterValue(model);
        r = treatmentProbability.getParameterValue(model);
        rho = samplingProbability.getParameterValue(model);
    }

    @Override
    public void updateLikelihoodModelValues(int model) {
        
        modelStartTime = modelStartTimes[model];

        if (model == 0) {
            previousP = 1.0;
        } else{
            previousP = p(model - 1, modelStartTime);
        }

        updateParameterValues(model);


        A = computeA(lambda, mu, psi);
        B = computeB(lambda, mu, psi, rho, A, previousP);
    }

    @Override
    public void updateGradientModelValues(int model) {

        modelStartTime = modelStartTimes[model];

        if (model == 0) {
            previousP = 1.0;
        } else{
            previousP = p(eAt_End);
        }

        updateParameterValues(model);

        A = computeA(lambda, mu, psi);
        B = computeB(lambda, mu, psi, rho, A, previousP);
        
        this.savedQ = Double.MIN_VALUE; // TODO What is this all about?
        this.partialQKnown = false;

        dACompute(dA);
        
        dBCompute(model, dB);

        if (numIntervals > 1 && model < numIntervals - 1) {

            double end = modelStartTimes[model + 1];
            double start = modelStartTimes[model];
            eAt_End = Math.exp(A * (end - start));
            dPCompute(model, end, start, eAt_End, dPModelEnd);
        }
        
        computedBCurrent = true;
    }

    @Override
    public double processInterval(int model, double tYoung, double tOld, int nLineages) {
        return nLineages * (logQ(model, tOld) - logQ(model, tYoung));
    }

    @Override
    public double processCoalescence(int model, double tOld) {
        return Math.log(lambda);
    }

    @Override
    public double processSampling(int model, double tOld) {

        double logSampProb;

        boolean sampleIsAtEventTime = tOld == modelStartTimes[model];
        boolean samplesTakenAtEventTime = rho > 0;

        if (sampleIsAtEventTime && samplesTakenAtEventTime) {
            logSampProb = Math.log(rho);
            if (model > 0) {
                logSampProb = Math.log(rho) + Math.log(r + (1.0 - r) * p(model,tOld));
            }
        } else {
            double logPsi = Math.log(psi);
            logSampProb = logPsi + Math.log(r + (1.0 - r) * p(model,tOld));
        }

        return logSampProb;
    }

    @Override
    public double processOrigin(int model, double rootAge) {
        if (originTime.getValue(0) < rootAge) {
            return Double.NaN;
        } else {
            double[] modelBreakPoints = getBreakPoints();
            double logL = 0.0;

            double origin = originTime.getValue(0);
            double intervalStart = rootAge;
            double segmentIntervalEnd = modelBreakPoints[model];

            while (origin >= segmentIntervalEnd) { // TODO Maybe it's >= ?
                logL += processModelSegmentBreakPoint(model, intervalStart, segmentIntervalEnd, 1);
                intervalStart = segmentIntervalEnd;
                ++model;
                updateLikelihoodModelValues(model);
                segmentIntervalEnd = modelBreakPoints[model];
            }
            if (intervalStart < origin) {
                logL += logQ(model, origin) - logQ(model, intervalStart);
            }

            return logL;
        }
    }

    public double calculateTreeLogLikelihood(Tree tree, Set<Taxon> exclude) {
        if (exclude.size() == 0) return calculateTreeLogLikelihood(tree);
        throw new RuntimeException("Not implemented!");
    }

    @Override
    public Citation.Category getCategory() {
        return Citation.Category.TREE_PRIORS;
    }

    @Override
    public String getDescription() {
        return "A (possibly episodic) serially-sampled birth-death model with the possibility of treatment (death when sampled) and intensive sampling.";
    }

    @Override
    public List<Citation> getCitations() {
        return Collections.singletonList(new Citation(
                new Author[]{
                        new Author("A", "Gavryushkina"),
                },
                "Bayesian Inference of Sampled Ancestor Trees",
                2014,
                "PLoS Computational Biology",
                10,
                1, 15,
                "10.1371/journal.pcbi.1003919"
        ));
    }

    private double g1(double eAt) {
        return  (1 + B) * eAt + (1 - B);
    }

    private double g2(double G1) {
        return A * (1 - 2 * (1 - B) / G1);
    }

    public double q(int model, double t) {
        double eAt = Math.exp(A * (t - modelStartTimes[model]));
        return q(eAt);
    }

    public double q(double eAt) {
        double sqrtDenominator = g1(eAt);
        return 4 * eAt / (sqrtDenominator * sqrtDenominator);
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

    private void dACompute(double[] dA) {
        dA[0] = (lambda - mu + psi) / A;
        dA[1] = (-lambda + mu + psi) / A;
        dA[2] = (lambda + mu + psi) / A;
    }

    private void dBCompute(int model, double[] dB) {

        for (int k = 0; k < model; ++k) {
            for (int p = 0; p < 4; p++) {
                dB[k * 4 + p] = -2 * (1 - rho) * lambda / A * dPModelEnd[k * 4 + p];
            }
        }

        double term1 = 1 - 2 * (1 - rho) * previousP;

        dB[model * 4 + 0] = (A * term1 - dA[0] * (term1 * lambda + mu + psi)) / (A * A);
        dB[model * 4 + 1] = (A - dA[1] * (term1 * lambda + mu + psi)) / (A * A);
        dB[model * 4 + 2] = (A - dA[2] * (term1 * lambda + mu + psi)) / (A * A);
        dB[model * 4 + 3] = 2 * lambda * previousP / A;
    }

    private void dPCompute(int model, double t, double intervalStart, double eAt, double[] dP) {

        double G1 = g1(eAt);

        double term1 = -A / lambda * ((1 - B) * (eAt - 1) + G1) / (G1 * G1);

        for (int k = 0; k  < model; k ++) {
            for (int p = 0; p < 4; p++) {
                dP[k * 4 + p] = term1 * dB[k * 4 + p];
            }
        }

        for (int p = 0; p < 3; ++p) { // TODO Only dG1[1] and dG2[2] are used
            double term2 = eAt * (1 + B) * dA[p] * (t - intervalStart) + (eAt - 1) * dB[model * 4 + p];
            dG2[p] = dA[p] - 2 * (G1 * (dA[p] * (1 - B) - dB[model * 4 + p] * A) - (1 - B) * term2 * A) / (G1 * G1);
        }

        double G2 = g2(G1);

        dP[model * 4 + 0] = (-mu - psi - lambda * dG2[0] + G2) / (2 * lambda * lambda);
        dP[model * 4 + 1] = (1 - dG2[1]) / (2 * lambda);
        dP[model * 4 + 2] = (1 - dG2[2]) / (2 * lambda);
        dP[model * 4 + 3] = -A / lambda * ((1 - B) * (eAt - 1) + G1) * dB[model * 4 + 3] / (G1 * G1);
    }

    private void dQCompute(int model, double t, double[] dQ) {

        double eAt = Math.exp(A * (t - modelStartTimes[model]));
        dQCompute(model, t, dQ, eAt);
    }

    private void dQCompute(int model, double t, double[] dQ, double eAt) {

        double dwell = t - modelStartTimes[model];
        double G1 = g1(eAt);

        double term1 = 8 * eAt;
        double term2 = G1 / 2 - eAt * (1 + B);
        double term3 = eAt - 1;
        double term4 = G1 * G1 * G1;
        double term5 = -term1 * term3 / term4;

        for (int k = 0; k < model; ++k) {
            for (int p = 0; p < 4; ++p) {
                dQ[k * 4 + p] = term5 * dB[k * 4 + p];
            }
        }

        double term6 = term1 / term4;
        double term7 = dwell * term2;

        dQ[model * 4 + 0] = term6 * (dA[0] * term7 - dB[model * 4 + 0] * term3);
        dQ[model * 4 + 1] = term6 * (dA[1] * term7 - dB[model * 4 + 1] * term3);
        dQ[model * 4 + 2] = term6 * (dA[2] * term7 - dB[model * 4 + 2] * term3);
        dQ[model * 4 + 3] = term5 * dB[model * 4 + 3];
    }

    @Override
    public void precomputeGradientConstants() {

        this.savedQ = Double.MIN_VALUE;
        this.partialQKnown = false;
    }

    public void processGradientModelSegmentBreakPoint(double[] gradient, int currentModelSegment,
                                                      double intervalStart, double intervalEnd, int nLineages) {
        double qStart;
        if (eAt_Old == 0) {
            qStart = q(currentModelSegment, intervalStart);
        } else {
            qStart = q(eAt_Old);
        }

        /*        double qStart = q(currentModelSegment, intervalStart);*/
        double qEnd = q(currentModelSegment, intervalEnd);

        dQCompute(currentModelSegment, intervalStart, dQStart);
        dQCompute(currentModelSegment, intervalEnd, dQEnd, eAt_End);

        for (int k = 0; k <= currentModelSegment; ++k) {
            for (int p = 0; p < 4; ++p) {
                gradient[genericIndex(k, p, numIntervals)] += nLineages *
                        (dQEnd[k * 4 + p] / qEnd - dQStart[k * 4 + p] / qStart);
            }
        }

        if ( samplingProbability.getValue(currentModelSegment + 1) > 0.0 && samplingProbability.getValue(currentModelSegment + 1) < 1.0) {
            // Add in probability of un-sampled lineages
            // We don't need this at t=0 because all lineages in the tree are sampled
            // TODO: check if we're right about how many lineages are actually alive at this time. Are we inadvertently over-counting or under-counting due to samples added at this _exact_ time?
            gradient[fractionIndex(currentModelSegment+1, numIntervals)] -= nLineages / (1.0 - samplingProbability.getValue(currentModelSegment + 1));
        }
    }

    @Override
    public void processGradientInterval(double[] gradient, int currentModelSegment,
                                        double intervalStart, double intervalEnd, int nLineages) {

        double[] partialQ_all_old = temp44;

        eAt_Old = Math.exp(A * (intervalEnd - modelStartTimes[currentModelSegment]));
        // dQCompute(currentModelSegment, intervalEnd, partialQ_all_old);
        dQCompute(currentModelSegment, intervalEnd, partialQ_all_old, eAt_Old);

        double[] partialQ_all_young;
        // double Q_Old = q(currentModelSegment, intervalEnd);
        double Q_Old = q(eAt_Old);
        double Q_young;
        if (this.savedQ != Double.MIN_VALUE) {
            Q_young = this.savedQ;
        } else {
            Q_young = q(currentModelSegment, intervalStart);
        }
        this.savedQ = Q_Old;

        if (partialQKnown) {
            // Only need 1 copy + 1 swap (instead of 2 copies)
            double[] tmp = partialQ;
            partialQ = temp33;
            temp33 = tmp;

            partialQ_all_young = temp33;
//            System.arraycopy(partialQ, 0, partialQ_all_young, 0, 4 * (currentModelSegment + 1));
        } else {
            partialQ_all_young = temp33;
            dQCompute(currentModelSegment, intervalStart, partialQ_all_young);
            partialQKnown = true;
        }
        System.arraycopy(partialQ_all_old, 0, partialQ, 0, 4 * (currentModelSegment + 1));

        for (int k = 0; k <= currentModelSegment; k++) {
            for (int p = 0; p < 4; p++) {
                gradient[k * 5 + p] += nLineages * (partialQ_all_old[k * 4 + p] / Q_Old - partialQ_all_young[k * 4 + p] / Q_young);
            }
        }
    }

    @Override
    public void processGradientSampling(double[] gradient, int currentModelSegment, double intervalEnd) {

        //boolean sampleIsAtPresent = intervalEnd <= 0;
        //boolean samplesTakenAtPresent = rho0 > 0;

        //TODO: need to confirm intensive sampling case is correct
        boolean sampleIsAtEventTime = intervalEnd == modelStartTimes[currentModelSegment];
        boolean samplesTakenAtEventTime = rho > 0;

        //if (sampleIsAtPresent && samplesTakenAtPresent) {
            //gradient[fractionIndex(0, numIntervals)] += 1 / rho;
        if (sampleIsAtEventTime && samplesTakenAtEventTime) {
            gradient[fractionIndex(currentModelSegment, numIntervals)] += 1 / rho; // TODO Need to test!
        } else {
            gradient[samplingIndex(currentModelSegment, numIntervals)] += 1 / psi;
        }

        if (!sampleIsAtEventTime | currentModelSegment > 0) {
            //double p_it = p(currentModelSegment, intervalEnd);
            if (intervalEnd == modelStartTimes[currentModelSegment]) {
                eAt_Old = Math.exp(A * (intervalEnd - modelStartTimes[currentModelSegment]));
            }

            double p_it = p(eAt_Old);

            gradient[treatmentIndex(currentModelSegment, numIntervals)] +=  (1 - p_it) / ((1 - r) * p_it + r);

            // double eAt = Math.exp(A * (intervalEnd - modelStartTimes[currentModelSegment]));

            dPCompute(currentModelSegment, intervalEnd, modelStartTimes[currentModelSegment], eAt_Old, this.dPIntervalEnd);

            double term1 = (1 - r) / ((1 - r) * p_it + r);

            for (int k = 0; k <= currentModelSegment; k++) {
                gradient[birthIndex(k, numIntervals)] += term1 * dPIntervalEnd[k * 4 + 0];
                gradient[deathIndex(k, numIntervals)] += term1 * dPIntervalEnd[k * 4 + 1];
                gradient[samplingIndex(k, numIntervals)] += term1 * dPIntervalEnd[k * 4 + 2];
                gradient[fractionIndex(k, numIntervals)] += term1 * dPIntervalEnd[k * 4 + 3];
            }
        }
    }

    @Override
    public void processGradientCoalescence(double[] gradient, int currentModelSegment, double intervalEnd) {
        gradient[birthIndex(currentModelSegment, numIntervals)] += 1 / lambda;
    }

    @Override
    public void processGradientOrigin(double[] gradient, int currentModelSegment, double totalDuration) {

        double origin = originTime.getValue(0);

        double[] modelBreakPoints = getBreakPoints();

        double intervalStart = totalDuration;

        while (origin >= modelBreakPoints[currentModelSegment]) { // TODO Maybe it's >= ?
            final double segmentIntervalEnd = modelBreakPoints[currentModelSegment];
            processGradientModelSegmentBreakPoint(gradient, currentModelSegment, intervalStart, segmentIntervalEnd, 1);
            intervalStart = segmentIntervalEnd;
            ++currentModelSegment;
            updateGradientModelValues(currentModelSegment);
        }
        if (intervalStart < origin) {
            double qOrigin = q(currentModelSegment, origin);
            double qIntervalStart = q(currentModelSegment, intervalStart);

            double[] dQOrigin = dQEnd;
            double[] dQIntervalStart = dQStart;

            dQCompute(currentModelSegment, origin, dQOrigin);
            dQCompute(currentModelSegment, intervalStart, dQIntervalStart);

            for (int k = 0; k <= currentModelSegment; k++) {
                gradient[birthIndex(k, numIntervals)] += dQOrigin[k * 4 + 0] / qOrigin - dQIntervalStart[k * 4 + 0] / qIntervalStart;
                gradient[deathIndex(k, numIntervals)] += dQOrigin[k * 4 + 1] / qOrigin - dQIntervalStart[k * 4 + 1] / qIntervalStart;
                gradient[samplingIndex(k, numIntervals)] += dQOrigin[k * 4 + 2] / qOrigin - dQIntervalStart[k * 4 + 2] / qIntervalStart;
                gradient[fractionIndex(k, numIntervals)] += dQOrigin[k * 4 + 3] / qOrigin - dQIntervalStart[k * 4 + 3] / qIntervalStart;
            }
        }


    }

    @Override
    public void logConditioningProbability(double[] gradient) {
        throw new RuntimeException("Cannot yet condition ESSBDP for gradient.");
    }

    @Override
    public int getGradientLength() { return 5 * numIntervals; }

    private static final boolean AOS = true;

    private static int genericIndex(final int n, final int p, final int stride) {
        if (AOS) {
            return n * 5 + p;
        } else {
            return p * stride + n;
        }
    }

    private static int birthIndex(final int n, final int stride) {
        if (AOS) {
            return n * 5;
        } else {
            return n;
        }
    }

    private static int deathIndex(final int n, final int stride) {
        if (AOS) {
            return n * 5 + 1;
        } else {
            return stride + n;
        }
    }

    private static int samplingIndex(final int n, final int stride) {
        if (AOS) {
            return n * 5 + 2;
        } else {
            return 2 * stride + n;
        }
    }

    private static int fractionIndex(final int n, final int stride) {
        if (AOS) {
            return n * 5 + 3;
        } else {
            return 3 * stride + n;
        }
    }

    private static int treatmentIndex(final int n, final int stride) {
        if (AOS) {
            return n * 5 + 4;
        } else {
            return 4 * stride + n;
        }
    }

    private static void swap(double[] array, final int i, final int j) {
        double tmp = array[i];
        array[i] = array[j];
        array[j] = tmp;
    }

    private static void transpose(double[] array, final int majorDim, final int minorDim) { // TODO untested
        int oldIndex = 0;
        for (int major = 0; major < majorDim; ++major) { // TODO Not cache-friendly, see https://stackoverflow.com/questions/5200338/a-cache-efficient-matrix-transpose-program
            for (int minor = major; minor < minorDim; ++minor) {
                final int newIndex = minor * minorDim + majorDim;
                swap(array, newIndex, oldIndex);
                ++oldIndex;
            }
        }
    }

    // Some useful JVM flags: -XX:+PrintCompilation -XX:+UnlockDiagnosticVMOptions -XX:+PrintInlining

    public abstract class ParameterPack implements Iterable<Parameter> {

        public class Constant extends ParameterPack {

            public Constant(Parameter birthRate,
                            Parameter deathRate,
                            Parameter samplingRate,
                            Parameter samplingProbability,
                            Parameter treatmentProbability,
                            Parameter originTime) {
                super(birthRate, deathRate, samplingRate, samplingProbability, treatmentProbability, originTime, null);
            }

            public double getBirthRate(int model) { return birthRate.getParameterValue(0); }

            public double getDeathRate(int model) { return deathRate.getParameterValue(0); }

            public double getSamplingRate(int model) { return samplingRate.getParameterValue(0); }

            public double getSamplingProbability(int model) { return samplingProbability.getParameterValue(0); }

            public double getTreatmentProbability(int model) { return treatmentProbability.getParameterValue(0); }
        }

        public class Episodic extends ParameterPack {

            public Episodic(Parameter birthRate,
                            Parameter deathRate,
                            Parameter samplingRate,
                            Parameter samplingProbability,
                            Parameter treatmentProbability,
                            Parameter originTime,
                            Parameter grid) {
                super(birthRate, deathRate, samplingRate, samplingProbability, treatmentProbability, originTime, grid);
            }

            public double getBirthRate(int model) { return birthRate.getParameterValue(model); }

            public double getDeathRate(int model) { return deathRate.getParameterValue(model); }

            public double getSamplingRate(int model) { return samplingRate.getParameterValue(model); }

            public double getSamplingProbability(int model) { return samplingProbability.getParameterValue(model); }

            public double getTreatmentProbability(int model) { return treatmentProbability.getParameterValue(model); }
        }

        final Parameter birthRate;
        final Parameter deathRate;
        final Parameter samplingRate;
        final Parameter samplingProbability;
        final Parameter treatmentProbability;
        final Parameter originTime;
        final Parameter grid;

        final int dim;

        final List<Parameter> parameterList = new ArrayList<>();

        private ParameterPack(Parameter birthRate,
                              Parameter deathRate,
                              Parameter samplingRate,
                              Parameter samplingProbability,
                              Parameter treatmentProbability,
                              Parameter originTime,
                              Parameter grid) {

            this.birthRate = birthRate;
            this.deathRate = deathRate;
            this.samplingRate = samplingRate;
            this.samplingProbability = samplingProbability;
            this.treatmentProbability = treatmentProbability;
            this.originTime = originTime;
            this.grid = grid;

            this.dim = grid == null ? 1 : grid.getDimension() + 1;

            addIfNotNull(birthRate, Double.POSITIVE_INFINITY, 0.0, dim);
            addIfNotNull(deathRate, Double.POSITIVE_INFINITY, 0.0, dim);
            addIfNotNull(samplingRate, Double.POSITIVE_INFINITY, 0.0, dim);
            addIfNotNull(samplingProbability, 1.0, 0.0, dim);
            addIfNotNull(treatmentProbability, 1.0, 0.0, dim);
            addIfNotNull(originTime, Double.POSITIVE_INFINITY, 0.0, 1);
        }

        public abstract double getBirthRate(int model);
        public abstract double getDeathRate(int model);
        public abstract double getSamplingRate(int model);
        public abstract double getSamplingProbability(int model);
        public abstract double getTreatmentProbability(int model);

        public double getOriginTime() {
            return originTime.getParameterValue(0);
        }

        public double[] getGrid() {
            double[] value = new double[dim + 1];
            value[0] = 0.0;
            for (int i = 1; i < dim; ++i) {
                value[i] = grid.getParameterValue(i - 1);
            }
            value[dim] = Double.POSITIVE_INFINITY;

            return value;
        }

        private void addIfNotNull(Parameter parameter, double upper, double lower, int dim) {
            if (parameter != null) {

                if (parameter.getDimension() != dim) {
                    throw new IllegalArgumentException("Parameter '" + parameter.getId() + "' has the wrong dimension");
                }

                parameterList.add(parameter);
                parameter.addBounds(new Parameter.DefaultBounds(upper, lower, dim));
            }
        }

        public boolean contains(Variable variable) {
            return parameterList.contains((Parameter) variable);
        }

        @Override
        public Iterator<Parameter> iterator() {
            return parameterList.iterator();
        }
    }
}
