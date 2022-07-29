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

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.util.Taxon;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;
import dr.util.Author;
import dr.util.Citable;
import dr.util.Citation;

import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Set;

/**
 * A phylogenetic birth-death-sampling model which includes serial sampling, sampling at present, and the possibility of treatmentProbability.
 */
public class NewBirthDeathSerialSamplingModel extends SpeciationModel implements SpeciationModelGradientProvider, Citable {
    // TODO should we pre-emptively put the combinatorial constant in here? It only really matters when inferring a tree where there may/may not be sampled ancestors

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
    double absTol = 1e-8;

    // there are numIntervals intervalStarts, implicitly intervalStarts[-1] == 0
    int numIntervals;
    double gridEnd;
    double[] intervalEnds = null;

    // TODO if we want to supplant other birth-death models, need an ENUM, and choice of options
    // Minimally, need survival of 1 lineage (passable default for SSBDP) and nTaxa (which is current option for non-serially-sampled BDP)
    private final boolean conditionOnSurvival;

    // useful constants we don't want to compute nTaxa times
    double A;
    double B;
    double pPrevious;
    double lambda;
    double mu;
    double psi;
    double r;
    double rho;

    private double[] savedGradient;
    private double savedQ;
    private double[] partialQ;
    private boolean partialQKnown;


    private double[] temp2;
    private double[] temp3;

    private double[] temp33;

    private double[] tempA;

    boolean computedBCurrent;
    private double[][] partialBCurrentPartialAll;

    private double[][] partialPPreviousPartialAll;

    private double[][] partialPCurrentPartialAll;



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

        this("NewBirthDeathSerialSamplingModel", birthRate, deathRate, serialSamplingRate, treatmentProbability, samplingProbability, originTime, condition, numIntervals, gridEnd, units);
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

        this.savedGradient = null;
//        this.savedTreeInterval = null;
        this.savedQ = Double.MIN_VALUE;
        this.partialQ = new double[4 * numIntervals];
        this.partialQKnown = false;

        this.tempA = new double[4];
        this.temp2 = new double[4];
        this.temp3 = new double[4];

        this.temp33 = new double[numIntervals * 4];

        this.gridEnd = gridEnd;
        this.numIntervals = numIntervals;
        setupTimeline();
    }

    // TODO should probably be replaced and brought in line with smoothSkygrid
    private void setupTimeline() {
        if (intervalEnds == null) {
            intervalEnds = new double[numIntervals];
        } else {
            Arrays.fill(intervalEnds, 0.0);
        }

        for (int idx = 0; idx <= numIntervals - 1 ; idx++) {
            intervalEnds[idx] = (idx+1) * (gridEnd / numIntervals);
        }
    }

    @Override
    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
//        fireModelChanged(variable);
    }

    /**
     * @param lambda   birth rate
     * @param mu   death rate
     * @param psi   proportion sampled at final time point
     * @param rho rate of sampling per lineage per unit time
     * @param t   time
     * @param eAt precomputed exp(A * (t - t_i))
     * @return the probability of no sampled descendants after time, t
     */
    // TODO make this take in the most recent break time as an argument or the model index so we can obtain said time
    // TODO do we really need 4 p functions?
    public static double p(int model, double lambda, double mu, double psi, double rho, double a, double b, double t, double ti, double eAt) {
        double eAt1B = eAt * (1.0 + b);
        return (lambda + mu + psi - a * ((eAt1B - (1.0 - b)) / (eAt1B + (1.0 - b)))) / (2.0 * lambda);
    }

    public static double p(int model, double lambda, double mu, double psi, double rho, double a, double b, double t, double ti) {
        double eAt = Math.exp(a * (t - ti));
        return p(model, lambda, mu, psi, rho, a, b, t, ti, eAt);
    }

    public double p(int model, double t, double eAt) {
        double ti = model == 0 ? 0 : intervalEnds[model-1];
        return p(model, lambda, mu, psi, rho, A, B, t, ti, eAt);
    }

    public double p(int model, double t) {
        double ti = model == 0 ? 0 : intervalEnds[model-1];
        return p(model, lambda, mu, psi, rho, A, B, t, ti);
    }

    /**
     * @param t   time
     * @return the probability of no sampled descendants after time, t
     */
    // TODO do we really need 4 logq functions?
    public double logq(int model, double t) {
        double ti = model == 0 ? 0 : intervalEnds[model-1];

        double At = A * (t - ti);
        double eA = Math.exp(At);
        double sqrtDenom = eA * (1.0 + B) + (1.0 - B);
        //System.out.println(At + Math.log(4.0) - 2.0 * Math.log(sqrtDenom));
        return At + Math.log(4.0) - 2.0 * Math.log(sqrtDenom);
    }
    public double logq(int model, double t, double eAt) {
        return Math.log(4.0 * eAt) - 2.0 * Math.log(eAt * (1 + B) + (1 - B));
    }

    // Named as per Gavryushkina et al 2014
    public static double computeA(double lambda, double mu, double psi) {
        return Math.abs(Math.sqrt(Math.pow(lambda - mu - psi, 2.0) + 4.0 * lambda * psi));
    }

    // Named as per Gavryushkina et al 2014
    public static double computeB(double lambda, double mu, double psi, double rho, double A, double pPrevious) {
        return ((1.0 - 2.0 * (1 - rho) * pPrevious) * lambda + mu + psi)/A;
    }

    @Override
    public double logConditioningProbability() {
        double logP = 0.0;
        if ( conditionOnSurvival ) {
            logP -= Math.log(1.0 - p(0, originTime.getParameterValue(0)));
        }
        return logP;
    }

    /**
     * Generic likelihood calculation
     *
     * @param tree the tree to calculate likelihood of
     * @return log-likelihood of density
     */
    public final double calculateTreeLogLikelihood(Tree tree) {
        // TODO deprecate this function, we only ever want to use the new loop
        updateModelValues();

        double logL = calculateUnconditionedTreeLogLikelihood(tree);

        double origin = originTime.getValue(0);
        if (origin < tree.getNodeHeight(tree.getRoot())) {
            return Double.NEGATIVE_INFINITY;
        }

        if ( conditionOnSurvival ) {
            logL -= Math.log(1.0 - p(0, origin));
        }

        return logL;
    }

    // Log-likelihood of tree without conditioning on anything
    public final double calculateUnconditionedTreeLogLikelihood(Tree tree) {

        // TODO deprecate this function, we only ever want to use the new loop

        double timeZeroTolerance = Double.MIN_VALUE;
        boolean noSamplingAtPresent = rho < Double.MIN_VALUE;

        double origin = originTime.getValue(0);
        if (origin < tree.getNodeHeight(tree.getRoot())) {
            return Double.NEGATIVE_INFINITY;
        }

        // extant leaves
        int n = 0;
        // extinct leaves
        int m = 0;
        // sampled ancestors
        int k = 0;

        for (int i = 0; i < tree.getExternalNodeCount(); i++) {
            NodeRef node = tree.getExternalNode(i);
            if (noSamplingAtPresent || tree.getNodeHeight(node) > timeZeroTolerance) {
                m += 1;
            } else {
                n += 1;
            }
        }

        if ((!noSamplingAtPresent) && n < 1) {
            throw new RuntimeException(
                    "Sampling fraction at time zero (rho) is >0 but there are no samples at time zero"
            );
        }

        double logL = 0.0;

        logL += (double)(n + m - 1) * Math.log(lambda);

        if ( k > 0 ) {
            logL += (double)(k) * Math.log(psi * (1.0 - r));
        }

        if (!noSamplingAtPresent) {
            logL += n * Math.log(rho);
        }

        logL -= logq(0, origin);

        for (int i = 0; i < tree.getInternalNodeCount(); i++) {
            double x = tree.getNodeHeight(tree.getInternalNode(i));
            logL -= logq(0, x);
        }

        double temp_eAt = 0;

        for (int i = 0; i < tree.getExternalNodeCount(); i++) {
            double y = tree.getNodeHeight(tree.getExternalNode(i));
            if (noSamplingAtPresent || y > timeZeroTolerance) {
                // TODO(change here)
                temp_eAt = Math.exp(A * y);
                logL += Math.log(psi * (r + (1.0 - r) * p(0, y, temp_eAt))) + logq(0, y, temp_eAt);
//                System.err.println("logq(y) = " + logq(y));
            }
        }

        return logL;
    }

    private double lambda(int i) {
        return birthRate.getValue(i);
    }

    private double mu(int i) {
        return deathRate.getValue(i);
    }

    private double psi(int i) {
        return serialSamplingRate.getValue(i);
    }

    private double rho(int i) {return samplingProbability.getValue(i); }

    public double r(int i) {
        return treatmentProbability.getValue(i);
    }

    @Override
    public double[] getBreakPoints() {
        return intervalEnds;
    }

    @Override
    public double processModelSegmentBreakPoint(int model, double intervalStart, double segmentIntervalEnd, int nLineages) {
        return nLineages * (logq(model, segmentIntervalEnd) - logq(model, intervalStart));
    }

    @Override
    public void updateModelValues(int model) {
        // TODO would it be excessive to also get the log-values for these? Would save O(n) calls to Math.log()
        this.savedQ = Double.MIN_VALUE;
        this.partialQKnown = false;
        if (model == 0) {
            pPrevious = 1.0;
        } else{
            pPrevious = p(model-1, intervalEnds[model-1]);
        }

        lambda = birthRate.getParameterValue(model);
        mu = deathRate.getParameterValue(model);
        psi = serialSamplingRate.getParameterValue(model);
        r = treatmentProbability.getParameterValue(model);
        rho = samplingProbability.getParameterValue(model);


        A = computeA(lambda, mu, psi);
        B = computeB(lambda, mu, psi, rho, A, pPrevious);

        partialApartialAll();

        computedBCurrent = false;
        // computedPPrevious = false;
        //partialBPreviousPartialAll = partialBCurrentPartialAll;

        if (partialBCurrentPartialAll == null) {
            partialPPreviousPartialAll = new double[numIntervals][4];
            partialPCurrentPartialAll = new double[numIntervals][4];
            partialBCurrentPartialAll = new double[numIntervals][4];
        }

        if (numIntervals > 1) {
            partialPPreviousPartialAll = partialPCurrentPartialAll;

            for (int k = 0; k <= model; k++) {
                double ti = model == 0 ? 0 : intervalEnds[model - 1];
                double eAt = Math.exp(A * (intervalEnds[model] - ti));
                partialPCurrentPartialAll[k] = partialPpartialAll(model, k, intervalEnds[model], eAt);
            }
        }
    }

    @Override
    public double processInterval(int model, double tYoung, double tOld, int nLineages) {
        return nLineages * (logq(model, tOld) - logq(model, tYoung));
    }

    @Override
    public double processCoalescence(int model, double tOld) {
//        return Math.log(lambda()); // TODO Notice the natural parameterization is `log lambda`
        return Math.log(lambda);
    }

    @Override
    public double processSampling(int model, double tOld) {

        double logSampProb = 0.0;


        boolean sampleIsAtPresent = tOld <= 0;
        boolean samplesTakenAtPresent = rho(0) >= 0;

        boolean sampleIsAtEventTime = Math.abs(tOld - intervalEnds[model]) <= 0;
        boolean samplesTakenAtEventTime = rho(model) >= 0;

        if (sampleIsAtPresent && samplesTakenAtPresent) {
            logSampProb = Math.log(rho);
        } else if (sampleIsAtEventTime && samplesTakenAtEventTime) {
            logSampProb = Math.log(rho);
        } else {
            double logPsi = Math.log(psi); // TODO Notice the natural parameterization is `log psi`
            logSampProb = logPsi + Math.log(r + (1.0 - r) * p(model,tOld));
        }

        return logSampProb;
    }

    @Override
    public double processOrigin(int model, double rootAge) {
        if (originTime.getValue(0) < rootAge) {
            return Double.NaN;
        } else {
            return (logq(model, originTime.getValue(0))) - logq(model, rootAge);
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

    // Material from `Gradient` class


    private double g1(int model, double t, double eAt) {
        return  (1 + B) * eAt + (1 - B);
    }

    private double g2(double t, double G1) {
        return A * (1 - 2 * (1 - B) / G1);
    }


/*    public double Q(double t){
        // TODO why the factor of 4 and inversion here?
        double eAt = Math.exp(A * t);
        return (2.0 * (1.0 - Math.pow(B,2.0)) + (1.0/eAt) * Math.pow((1.0 - B),2.0) + eAt * Math.pow(1.0 + B,2.0));
    }*/

    //  todo: avoid calculating g1 in q
    public double q(int model, double t) {
        double ti  = model == 0 ? 0 : intervalEnds[model-1];;

        double eA = Math.exp(A * (t - ti));
        double sqrtDenom = eA * (1.0 + B) + (1.0 - B);

        return (4.0 * eA)/(Math.pow(sqrtDenom,2.0));
    }


    private void partialApartialLambda(int idx1, double lambda, double mu, double psi) {
        this.tempA[idx1] = (lambda - mu + psi) / A;
    }

    private void partialApartialMu(int idx1, double lambda, double mu, double psi) {
        this.tempA[idx1] = (-lambda + mu + psi) / A;
    }

    private void partialApartialPsi(int idx1, double lambda, double mu, double psi) {
        this.tempA[idx1] = (lambda + mu + psi) / A;
    }

    private void partialApartialRho(int idx1) {
        this.tempA[idx1] = 0;
    }

    private void partialApartialAll() {
        partialApartialLambda(0,  lambda, mu, psi);
        partialApartialMu(1, lambda, mu, psi);
        partialApartialPsi(2, lambda, mu, psi);
        partialApartialRho(3);
    }

    // Gradient w.r.t. Rho
/*    private void partialC1C2partialRho(int idx1, int idx2) {
        // double lambda = lambda();

//        double[] partialC1C2 = new double[2];
        this.temp1[idx1] = 0;
        this.temp1[idx2] = 2 * lambda / A;

//        return partialC1C2;
    }

    private void partialC1C2partialMu(int idx1, int idx2, double lambda, double mu, double psi, double rho) {
        // c1 == constants[0], c2 == constants[1]
//        double[] constants = getConstants();
        // double lambda = lambda();
        // double mu = mu();
        // double psi = psi();
        // double rho = rho();

//        double[] partialC1C2 = new double[2];
        this.temp1[idx1] = (-lambda + mu + psi) / A;
        this.temp1[idx2] = (A + (lambda - mu - 2 * lambda * rho - psi) * this.temp1[idx1]) / (A * A);

//        return partialC1C2;
    }

    private void partialC1C2partialLambda(int idx1, int idx2, double lambda, double mu, double psi, double rho) {
        // c1 == constants[0], c2 == constants[1]
//        double[] constants = getConstants();
        // double lambda = lambda();
        // double mu = mu();
        // double psi = psi();
        // double rho = rho();

//        double[] partialC1C2 = new double[2];
        this.temp1[idx1] = (lambda - mu + psi) / A;
        this.temp1[idx2] = ((2*rho - 1)*A - (-lambda + mu + 2 * lambda * rho + psi) * this.temp1[idx1]) / (A * A);

//        return partialC1C2;
    }

    private void partialC1C2partialPsi(int idx1, int idx2, double lambda, double mu, double psi, double rho) {
        // c1 == constants[0], c2 == constants[1]
//        double[] constants = getConstants();
        // double lambda = lambda();
        // double mu = mu();
        // double psi = psi();
        // double rho = rho();

//        double[] partialC1C2 = new double[2];
        this.temp1[idx1] = (lambda + mu + psi) / A;
        this.temp1[idx2] = (A + (lambda - mu - 2 * lambda * rho - psi) * this.temp1[idx1]) / (A * A);

//        return partialC1C2;
    }*/

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


/*    // gradients for all
    // (lambda, mu, psi, rho)
    public void partialC1C2partialAll() {
        // double lambda = lambda();
        // double mu = mu();
        // double psi = psi();
        // double rho = rho();

        partialC1C2partialLambda(0, 1, lambda, mu, psi, rho);
        partialC1C2partialMu(2, 3, lambda, mu, psi, rho);
        partialC1C2partialPsi(4, 5, lambda, mu, psi, rho);
        partialC1C2partialRho(6, 7);
//        partialC1C2_all[0] = partialC1C2partialLambda();
//        partialC1C2_all[1] = partialC1C2partialMu();
//        partialC1C2_all[2] = partialC1C2partialPsi();
//        partialC1C2_all[3] = partialC1C2partialRho();
//        return this.temp1;
    }

    // (lambda, mu, psi, rho)
//    public double[] partialQpartialAll(double t) {
//        double[] buffer = new double[4];
//        return partialQpartialAll(buffer, t);
//    }*/

    private double[] partialBpartialAll(int i, int k, double lambda, double mu, double psi, double rho) {
        if (computedBCurrent) {
            return partialBCurrentPartialAll[k];
        }
        double[] partialB = new double[4];
        if (k == i) {
            double[] partialA = tempA;
            double temp = 1 - 2 * (1 - rho) * pPrevious;
            partialB[0] = (A * temp - partialA[0] * (temp * lambda + mu + psi)) / (A* A);
            partialB[1] = (A - partialA[1] * (temp * lambda + mu+ psi)) / (A * A);
            partialB[2] = (A - partialA[2] * (temp * lambda + mu + psi)) / (A * A);
            partialB[3] = 2 * lambda * pPrevious / A;
        } else if (k < i) {
            for(int n = 0; n < 4; n++) {
                partialB[n] = -2 * (1 - rho) * lambda / A * partialPPreviousPartialAll[k][n];
            }
        }
        // else: k > i, all zero
        partialBCurrentPartialAll[k] = partialB;
        if (k == i) {
            computedBCurrent = true;
        }
        return partialB;
    }

/*
    public double[] partialQpartialAll(double[] partialQ_all, double t) {
//        double[] constants = getConstants();

        double eAt = Math.exp(-A * t);

//        double v = Math.exp(A * t) * (1 + B) - eAt * (1 - B) - 2 * B;
        double v = (1 + B) / eAt - eAt * (1 - B) - 2 * B;
        double v1 = (1 + B) /eAt * (1 + B) - eAt * (1 - B) * (1 - B);

        partialC1C2partialAll();

        double[] partialC1C2_all = this.temp1;

//        double[] partialQ_all = new double[4];
       // Arrays.fill(partialQ_all, 0.0);
        for (int i = 0; i < 4; ++i) {
            partialQ_all[i] = t * partialC1C2_all[i*2+0] * v1;
            partialQ_all[i] += 2 * partialC1C2_all[i*2+1] * v;
        }
        return partialQ_all;
    }
*/

    // (lambda, mu, psi, rho)
    public double[] partialG2partialAll(int model, double t, double eAt, double G1, double[] partialBStored) {

//        double eAt = Math.exp(A * (t-ti));
        double ti = model == 0 ? 0 : intervalEnds[model-1];
        double[] partialA = tempA;
        double[] partialB = partialBStored;
        double[] partialG2_all = temp2; // new double[4];
        double[] temp3 = new double[3];

        for (int n = 0; n < 3; ++n) {
            temp3[n] = eAt * (1 + B) * partialA[n] * (t - ti) + (eAt - 1) * partialB[n];
            double partialG2 = partialA[n] - 2 * (G1 * (partialA[n] * (1 - B) - partialB[n] * A) - (1 - B) * temp3[n] * A) / (G1 * G1);
            partialG2_all[n] = partialG2;
            }

        partialG2_all[3] = 0; // w.r.t. rho
        return partialG2_all;
    }

    // (lambda, mu, psi, rho)
    public double[] partialPpartialAll(int model, int k, double t, double eAt) {

        double[] partialP_all = new double[4];
        double[] partialB = partialBpartialAll(model, k, lambda, mu, psi, rho);
        // double G1 = g1(t);
        double G1 = g1(model, t, eAt);
        if (k == model) {
            double[] partialG2_all = partialG2partialAll(model, t, eAt, G1, partialB);
            double G2 = g2(t, G1);

                // lambda
            partialP_all[0] = (-mu - psi - lambda * partialG2_all[0] + G2) / (2 * lambda * lambda);
                // mu
            partialP_all[1] = (1 - partialG2_all[1]) / (2 * lambda);
                // psi
            partialP_all[2] = (1 - partialG2_all[2]) / (2 * lambda);
                // rho
            partialP_all[3] = -A / lambda * ((1 - B) * (eAt - 1) + G1) * partialB[3] / Math.pow(G1, 2);

        }  else if (k < model) {
            for (int n = 0; n < 4; n++) {
                partialP_all[n] = -A / lambda * ((1 - B) * (eAt - 1) + G1 ) * partialB[n] / Math.pow(G1,2);
            }
        }

        return partialP_all;
    }

    private double[] partialqpartialAll(int model, int k, double t, double[] partialq_all) {
        if (model < k) { return partialq_all;}
        double[] partialA = tempA;
        double tempA1;
        double ti = model == 0 ? 0 : intervalEnds[model-1];
        double eAt = Math.exp(A * (t - ti));
        double G1 = g1(model, t, eAt);
        double[] partialB = partialBpartialAll(model, k, lambda, mu, psi, rho);
        for (int n = 0; n < 3; n++) {
            if (model == k) {
                tempA1 = (t - ti) * partialA[n] * (G1 / 2 - eAt * (1 + B));
                //tempA2 = (t - ti) * partialA[n] * temp_exp * (1 + Bi[i]);
            } else {
                tempA1 = 0;
                //tempA2 = 0;
            }
            partialq_all[n] = 8 * eAt * (tempA1 - partialB[n] * (eAt - 1)) / (Math.pow(G1, 3));//(tempA1 - 8 * temp_exp * temp2 * (tempA2 + partialB[n] * (temp_exp - 1))) / (Math.pow(temp2, 4));
        }
        partialq_all[3] = -8 * eAt * partialB[3] * (eAt - 1) / (Math.pow(G1, 3));
        return partialq_all;
    }
    @Override
    public void precomputeGradientConstants() {
        this.savedQ = Double.MIN_VALUE;
//        this.savedPartialQ = null;
        this.partialQKnown = false;
        partialPPreviousPartialAll = new double[numIntervals][4];
        partialPCurrentPartialAll = new double[numIntervals][4];
        partialBCurrentPartialAll = new double[numIntervals][4];
    }

    @Override
    public void processGradientModelSegmentBreakPoint(double[] gradient, int currentModelSegment, double intervalStart, double segmentIntervalEnd, int nLineages) {
        double[] partialqStart = new double[4];
        double[] partialqEnd = new double[4];
        double qStart;
        double qEnd;
        for (int k = 0; k <= currentModelSegment; k++) {
            partialqpartialAll(currentModelSegment, k, intervalStart, partialqStart);
            partialqpartialAll(currentModelSegment, k, segmentIntervalEnd, partialqEnd);
            qStart = q(currentModelSegment, intervalStart);
            qEnd = q(currentModelSegment, segmentIntervalEnd);

            for (int n = 0; n < 4; n++) {
                // assume lambda_1, lambda_2, ..., mu_1, mu_2, ...
                gradient[k * 5 + n] += nLineages * (partialqEnd[n] / qEnd - partialqStart[n] / qStart);
                //gradient[n * numIntervals + k] += nLineages * (partialqEnd[n] / qEnd - partialqStart[n] / qStart);
            }
        }
    }

    @Override
    public void processGradientInterval(double[] gradient, int currentModelSegment, double intervalStart, double intervalEnd, int nLineages) {
        double tOld = intervalEnd;
        double tYoung = intervalStart;
        double[] partialQ_all_old = new double[numIntervals*4];
        double[] temp = new double[4];
        for (int k = 0; k <= currentModelSegment; k++) {
            // temp = partialqpartialAll(currentModelSegment, k, tOld, temp2);
            partialqpartialAll(currentModelSegment, k, tOld, temp);
            System.arraycopy(temp, 0, partialQ_all_old, 4*k, 4);
        }
        double[] partialQ_all_young;
        double Q_Old = q(currentModelSegment, tOld);
        double Q_young;
        if (this.savedQ != Double.MIN_VALUE) {
            Q_young = this.savedQ;
        } else {
            Q_young = q(currentModelSegment, tYoung);
        }
        this.savedQ = Q_Old;

        partialQ_all_young = temp33;
        if (partialQKnown) {
            // TODO: Maybe change numIntervals to currentModelSegment
            System.arraycopy(partialQ, 0, partialQ_all_young, 0, 4*numIntervals);
        } else {
            for (int k = 0; k <= currentModelSegment; k++) {
                temp = partialqpartialAll(currentModelSegment, k, tYoung, temp3);
                System.arraycopy(temp, 0, partialQ_all_young, 4*k, 4);
            }
            //System.arraycopy(partialQ_all_young, 0, savedPartialQ, 0, 4);
            partialQKnown = true;
        }
        // TODO: Maybe change numIntervals to currentModelSegment
        System.arraycopy(partialQ_all_old, 0, partialQ, 0, 4*numIntervals);

//        if (this.savedPartialQ != null) {
//            partialQ_all_young = this.savedPartialQ;
//        }
//        else {
//            partialQ_all_young = partialQpartialAll(tYoung);
//        }
//        this.savedPartialQ = partialQ_all_old;
        for (int k = 0; k <= currentModelSegment; k++) {
            for (int n = 0; n < 4; n++) {
                gradient[k * 5 + n] += nLineages * (partialQ_all_old[k*4+ n] / Q_Old - partialQ_all_young[k*4 + n] / Q_young);
                //gradient[n * numIntervals + k] += nLineages * (partialqEnd[n] / qEnd - partialqStart[n] / qStart);
            }
        }
    }

    @Override
    public void processGradientSampling(double[] gradient, int currentModelSegment, double intervalEnd) {

        boolean sampleIsAtPresent = intervalEnd <= 0;
        boolean samplesTakenAtPresent = rho(0) >= 0;

        boolean sampleIsAtEventTime = Math.abs(intervalEnd - intervalEnds[currentModelSegment]) <= 0;
        boolean samplesTakenAtEventTime = rho(currentModelSegment) >= 0;

        if (sampleIsAtPresent && samplesTakenAtPresent) {
            gradient[3] += 1 / rho;
            //gradient[3*numIntervals] += 1 / rho(0);
        } else if (sampleIsAtEventTime && samplesTakenAtEventTime) {
            gradient[3 + 5*currentModelSegment] += 1 / rho;
            //gradient[3*numIntervals + currentModelSegment] += 1 / rho(currentModelSegment);
        } else {
            // double logPsi = Math.log(serialSamplingRate(currentModelSegment)); // TODO Notice the natural parameterization is `log psi`
            gradient[2 + 5*currentModelSegment] += 1 / psi;
            //gradient[2*numIntervals + currentModelSegment] += 1 / psi(currentModelSegment);
            double p_it = p(currentModelSegment, intervalEnd);
            gradient[4 + 5*currentModelSegment] += (1 - p_it) / ((1-r)*p_it + r);
            //gradient[4*numIntervals + currentModelSegment] += (1 - p_it) / ((1-r)*p_it + r);
            double[] partialP;
            for(int k = 0; k <= currentModelSegment; k++) {
                double ti = currentModelSegment == 0 ? 0 : intervalEnds[currentModelSegment-1];
                double eAt = Math.exp(A * (intervalEnd - ti));
                partialP = partialPpartialAll(currentModelSegment, k, intervalEnd,  eAt);
                for(int n = 0; n < 4; n++) {
                    gradient[n + 5*k] += (1 - r) / ((1 - r) * p_it + r) * partialP[n];
                    //gradient[n*numIntervals + k] += (1 - r) / ((1 - r) * p_it + r) * partialP[n];
                }
            }
        }

    }

    @Override
    public void processGradientCoalescence(double[] gradient, int currentModelSegment, double intervalEnd) {
        gradient[currentModelSegment*5] += 1 / lambda;
    }

    @Override
    public void processGradientOrigin(double[] gradient, int currentModelSegment, double totalDuration) {
        double origin = originTime.getValue(0);
        double[] partialqOrigin;
        double[] partialqRoot;
        double qOrigin;
        double qRoot;

        for (int k = 0; k <= currentModelSegment; k++) {
            partialqOrigin = partialqpartialAll(currentModelSegment, k, origin, temp2);
            partialqRoot = partialqpartialAll(currentModelSegment, k, totalDuration, temp3);
            qOrigin = q(currentModelSegment, origin);
            qRoot = q(currentModelSegment, totalDuration);

            for (int n = 0; n < 4; n++) {
                gradient[k * 5 + n] += (partialqOrigin[n] / qOrigin - partialqRoot[n] / qRoot);
                //gradient[n * numIntervals + k] += (partialqOrigin[n] / qOrigin - partialqRoot[n] / qRoot);
            }
        }

    }

    @Override
    public void logConditioningProbability(double[] gradient) {
    }

    @Override
    public int getGradientLength() { return 5 * numIntervals; }

    private static final boolean AOS = true;

    private static final int birthIndex(final int n, final int stride) {
        if (AOS) {
            return n * stride + 0;
        } else {
            return 0 * stride + n;
        }
    }

    private static final int deathIndex(final int n, final int stride) {
        if (AOS) {
            return n * stride + 1;
        } else {
            return 1 * stride + n;
        }
    }

    private static final int samplingIndex(final int n, final int stride) {
        if (AOS) {
            return n * stride + 2;
        } else {
            return 2 * stride + n;
        }
    }

    private static final int fractionIndex(final int n, final int stride) {
        if (AOS) {
            return n * stride + 3;
        } else {
            return 3 * stride + n;
        }
    }

    private static final int treatmentIndex(final int n, final int stride) {
        if (AOS) {
            return n * stride + 4;
        } else {
            return 4 * stride + n;
        }
    }



    private static void swap(double[] array, final int i, final int j) {
        double tmp = array[i];
        array[i] = array[j];
        array[j] = tmp;
    }

    private static final void transpose(double[] array, final int majorDim, final int minorDim) { // TODO untested
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
}
