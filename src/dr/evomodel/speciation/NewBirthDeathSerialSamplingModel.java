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

    // extant sampling proportion
    Parameter samplingFractionAtPresent;

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

    // TODO if we want to supplant other birth-death models, need an ENUM, and choice of options
    // Minimally, need survival of 1 lineage (passable default for SSBDP) and nTaxa (which is current option for non-serially-sampled BDP)
    private final boolean conditionOnSurvival;

    // useful constants we don't want to compute nTaxa times
    private double C1 = Double.NEGATIVE_INFINITY;
    private double C2 = Double.NEGATIVE_INFINITY;

    private double[] savedGradient;
    private double savedQ;
    private double[] partialQ;
    private boolean partialQKnown;

    private double[][] temp1;
    private double[] temp2;
    private double[] temp3;



    public NewBirthDeathSerialSamplingModel(
            Parameter birthRate,
            Parameter deathRate,
            Parameter serialSamplingRate,
            Parameter treatmentProbability,
            Parameter samplingFractionAtPresent,
            Parameter originTime,
            boolean condition,
            Type units) {

        this("NewBirthDeathSerialSamplingModel", birthRate, deathRate, serialSamplingRate, treatmentProbability, samplingFractionAtPresent, originTime, condition, units);
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
            Parameter samplingFractionAtPresent,
            Parameter originTime,
            boolean condition,
            Type units) {

        super(modelName, units);

        this.birthRate = birthRate;
        addVariable(birthRate);
        birthRate.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, 1));

        this.deathRate = deathRate;
        addVariable(deathRate);
        deathRate.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, 1));

        this.serialSamplingRate = serialSamplingRate;
        addVariable(serialSamplingRate);
        serialSamplingRate.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, 1));

        this.samplingFractionAtPresent = samplingFractionAtPresent;
        addVariable(samplingFractionAtPresent);
        samplingFractionAtPresent.addBounds(new Parameter.DefaultBounds(1.0, 0.0, 1));

        this.treatmentProbability = treatmentProbability;
        addVariable(treatmentProbability);
        treatmentProbability.addBounds(new Parameter.DefaultBounds(1.0, 0.0, 1));

        this.conditionOnSurvival = condition;

        this.originTime = originTime;
        if (originTime != null) {
            addVariable(originTime);
            originTime.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, 1));
        }

        this.savedGradient = null;
//        this.savedTreeInterval = null;
        this.savedQ = Double.MIN_VALUE;
        this.partialQ = new double[4];
        this.partialQKnown = false;

        this.temp1 = new double[4][2];
        this.temp2 = new double[4];
        this.temp3 = new double[4];
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
     * @return the probability of no sampled descendants after time, t
     */
    public static double p0(double lambda, double mu, double psi, double rho, double c1, double c2, double t) {

        double expc1trc2 = Math.exp(-c1 * t) * (1.0 - c2);

        // Stadler 2011 p 349
        return (lambda + mu + psi + c1 * ((expc1trc2 - (1.0 + c2)) / (expc1trc2 + (1.0 + c2)))) / (2.0 * lambda);
    }

    /**
     * @param t   time
     * @return the probability of no sampled descendants after time, t
     */
    public static double logq(double c1, double c2, double t) {
        double expC1t = Math.exp(c1 * t);
        double q = 4.0/(2.0 * (1.0 - Math.pow(c2,2.0)) + (1.0/expC1t) * Math.pow((1.0 - c2),2.0) + expC1t * Math.pow(1.0 + c2,2.0));
        return Math.log(q);
    }

    private static double c1(double lambda, double mu, double psi) {
        return Math.abs(Math.sqrt(Math.pow(lambda - mu - psi, 2.0) + 4.0 * lambda * psi));
    }

    private static double c2(double lambda, double mu, double psi, double rho) {
        return -(lambda - 2.0 * rho * lambda - mu - psi)/c1(lambda, mu, psi);
    }

    public void precomputeConstants() {
        this.C1 = c1(lambda(), mu(), psi());
        this.C2 = c2(lambda(), mu(), psi(), rho());
    }

    public double getC1() { return C1; }

    public double getC2() { return C2; }

//    public double[] getConstants() {
//        double[] constants = {C1, C2};
//        return constants;
//    }

    public double p0(double t) {
        return p0(lambda(), mu(), psi(), rho(), C1, C2, t);
    }

    public double logq(double t) {
        return logq(C1, C2, t);
    }

    double lambda() {
        return birthRate.getParameterValue(0);
    }

    double mu() {
        return deathRate.getParameterValue(0);
    }

    double psi() {
        return serialSamplingRate.getParameterValue(0);
    }

    double r() {
        return treatmentProbability.getParameterValue(0);
    }

    double rho() {
        return samplingFractionAtPresent.getParameterValue(0);
    }

    @Override
    public double logConditioningProbability() {
        double logP = 0.0;
        if ( conditionOnSurvival ) {
            logP -= Math.log(1.0 - p0(originTime.getParameterValue(0)));
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
        precomputeConstants();

        double logL = calculateUnconditionedTreeLogLikelihood(tree);

        double origin = originTime.getValue(0);
        if (origin < tree.getNodeHeight(tree.getRoot())) {
            return Double.NEGATIVE_INFINITY;
        }

        if ( conditionOnSurvival ) {
            logL -= Math.log(1.0 - p0(origin));
        }

        return logL;
    }

    // Log-likelihood of tree without conditioning on anything
    public final double calculateUnconditionedTreeLogLikelihood(Tree tree) {

        double lambda = lambda();
        double mu = mu();
        double psi = psi();
        double r = r();
        double rho = rho();

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

        logL -= logq(origin);

        for (int i = 0; i < tree.getInternalNodeCount(); i++) {
            double x = tree.getNodeHeight(tree.getInternalNode(i));
            logL -= logq(x);
        }

        for (int i = 0; i < tree.getExternalNodeCount(); i++) {
            double y = tree.getNodeHeight(tree.getExternalNode(i));
            if (noSamplingAtPresent || y > timeZeroTolerance) {
                logL += Math.log(psi * (r + (1.0 - r) * p0(y))) + logq(y);
//                System.err.println("logq(y) = " + logq(y));
            }
        }

        return logL;
    }

    @Override
    public double[] getBreakPoints() {
        return new double[] { Double.POSITIVE_INFINITY };
    }

    @Override
    public double processInterval(int model, double tYoung, double tOld, int nLineages) {
        return nLineages * (logq(tOld) - logq(tYoung));
    }

    @Override
    public double processOrigin(int model, double rootAge) {
        return (logq(originTime.getValue(0))) - logq(rootAge);
    }

    @Override
    public double processCoalescence(int model, double tOld) {
        return Math.log(lambda()); // TODO Notice the natural parameterization is `log lambda`
    }

    @Override
    public double processSampling(int model, double tOld) {

        double logPsi = Math.log(psi()); // TODO Notice the natural parameterization is `log psi`
        double r = r();
        double logRho = Math.log(rho()); // TODO Notice the natural parameterization is `log rho`

        double timeZeroTolerance = Double.MIN_VALUE;
        boolean noSamplingAtPresent = rho() < Double.MIN_VALUE;

        if (noSamplingAtPresent || tOld > timeZeroTolerance) {
            return logPsi + Math.log(r + (1.0 - r) * p0(tOld));
        } else {
            return logRho;
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
        return "A serially-sampled birth-death model with the possibility of treatment and sampling at present.";
    }

    @Override
    public List<Citation> getCitations() {
        return Collections.singletonList(new Citation(
                new Author[]{
                        new Author("T", "Gernhard"),
                },
                "The conditioned reconstructed process",
                2008,
                "Journal of Theoretical Biology",
                253,
                769, 778,
                "10.1016/j.jtbi.2008.04.005"
        ));
    }

    // Material from `Gradient` class

    private double g1(double t) {
//        double[] constants = getConstants();
        double C1 = getC1();
        double C2 = getC2();
        double G1 = Math.exp(-C1 * t) * (1 - C2) + (1 + C2);
        return G1;
    }

    private double g2(double t) {
//        double[] constants = getConstants();
        double C1 = getC1();
        double C2 = getC2();
        double G1 = g1(t);
        double G2 = C1 * (1 - 2 * (1 + C2) / G1);
        return G2;
    }

    private double Q(double t){
        return 4*(1/Math.exp(logq(t)));
    }

    // Gradient w.r.t. Rho
    private void partialC1C2partialRho(double[] partialC1C2) {
        // c1 == constants[0], c2 == constants[1]
//        double[] constants = getConstants();
        double lambda = lambda();
        double C1 = getC1();

//        double[] partialC1C2 = new double[2];
        partialC1C2[0] = 0;
        partialC1C2[1] = 2 * lambda / C1;

//        return partialC1C2;
    }

    private void partialC1C2partialMu(double[] partialC1C2) {
        // c1 == constants[0], c2 == constants[1]
//        double[] constants = getConstants();
        double lambda = lambda();
        double mu = mu();
        double psi = psi();
        double rho = rho();
        double C1 = getC1();

//        double[] partialC1C2 = new double[2];
        partialC1C2[0] = (-lambda + mu + psi) / C1;
        partialC1C2[1] = (C1 + (lambda - mu - 2 * lambda * rho - psi) * partialC1C2[0]) / (C1 * C1);

//        return partialC1C2;
    }

    private void partialC1C2partialLambda(double[] partialC1C2) {
        // c1 == constants[0], c2 == constants[1]
//        double[] constants = getConstants();
        double lambda = lambda();
        double mu = mu();
        double psi = psi();
        double rho = rho();
        double C1 = getC1();

//        double[] partialC1C2 = new double[2];
        partialC1C2[0] = (lambda - mu + psi) / C1;
        partialC1C2[1] = ((2*rho - 1)*C1 - (-lambda + mu + 2 * lambda * rho + psi) * partialC1C2[0]) / (C1 * C1);

//        return partialC1C2;
    }

    private void partialC1C2partialPsi(double[] partialC1C2) {
        // c1 == constants[0], c2 == constants[1]
//        double[] constants = getConstants();
        double lambda = lambda();
        double mu = mu();
        double psi = psi();
        double rho = rho();
        double C1 = getC1();

//        double[] partialC1C2 = new double[2];
        partialC1C2[0] = (lambda + mu + psi) / C1;
        partialC1C2[1] = (C1 + (lambda - mu - 2 * lambda * rho - psi) * partialC1C2[0]) / (C1 * C1);

//        return partialC1C2;
    }

    @Override
    public Parameter getSamplingProbabilityParameter() {
        return samplingFractionAtPresent;
    }

    @Override
    public double[] getSamplingProbabilityGradient(Tree tree, NodeRef node) {
        double[] result = new double[1];
        result[0] = getAllGradient(tree, node)[3];
        return result;
    }

    @Override
    public Parameter getDeathRateParameter() {
        return deathRate;
    }

    @Override
    public double[] getDeathRateGradient(Tree tree, NodeRef node) {
        double[] result = new double[1];
        result[0] = getAllGradient(tree, node)[1];
        return result;
    }


    @Override
    public Parameter getBirthRateParameter() {
        return birthRate;
    }

    @Override
    public double[] getBirthRateGradient(Tree tree, NodeRef node) {
        double[] result = new double[1];
        result[0] = getAllGradient(tree, node)[0];
        return result;
    }

    @Override
    public Parameter getSamplingRateParameter() {
        return serialSamplingRate;
    }

    @Override
    public double[] getSamplingRateGradient(Tree tree, NodeRef node) {
        double[] result = new double[1];
        result[0] = getAllGradient(tree, node)[2];
        return result;
    }


    @Override
    public Parameter getTreatmentProbabilityParameter() {
        return treatmentProbability;
    }

    @Override
    public double[] getTreatmentProbabilityGradient(Tree tree, NodeRef node) {
        double[] result = new double[1];
        result[0] = getAllGradient(tree, node)[4];
        return result;
    }

    // gradients for all
    // (lambda, mu, psi, rho)
    public double[][] partialC1C2partialAll(double[][] partialC1C2_all) {
        partialC1C2partialLambda(partialC1C2_all[0]);
        partialC1C2partialMu(partialC1C2_all[1]);
        partialC1C2partialPsi(partialC1C2_all[2]);
        partialC1C2partialRho(partialC1C2_all[3]);
//        partialC1C2_all[0] = partialC1C2partialLambda();
//        partialC1C2_all[1] = partialC1C2partialMu();
//        partialC1C2_all[2] = partialC1C2partialPsi();
//        partialC1C2_all[3] = partialC1C2partialRho();
        return partialC1C2_all;
    }

    // (lambda, mu, psi, rho)
//    public double[] partialQpartialAll(double t) {
//        double[] buffer = new double[4];
//        return partialQpartialAll(buffer, t);
//    }

    public double[] partialQpartialAll(double[] partialQ_all, double t) {
//        double[] constants = getConstants();
        double C1 = getC1();
        double C2 = getC2();

        double expC1t = Math.exp(-C1 * t);

//        double v = Math.exp(C1 * t) * (1 + C2) - expC1t * (1 - C2) - 2 * C2;
        double v = (1 + C2) / expC1t - expC1t * (1 - C2) - 2 * C2;
        double v1 = (1 + C2) /expC1t * (1 + C2) - expC1t * (1 - C2) * (1 - C2);

        double[][] partialC1C2_all = partialC1C2partialAll(temp1);

//        double[] partialQ_all = new double[4];
        Arrays.fill(partialQ_all, 0.0);
        for (int i = 0; i < 4; ++i) {
            partialQ_all[i] += t * partialC1C2_all[i][0] * v1;
            partialQ_all[i] += 2 * partialC1C2_all[i][1] * v;
        }
        return partialQ_all;
    }

    // (lambda, mu, psi, rho)
    public double[] partialG2partialAll(double t, double expC1t) {
//        double[] constants = getConstants();
        double C1 = getC1();
        double C2 = getC2();

//        double expC1t = Math.exp(-C1 * t);

        double[][] partialC1C2_all = partialC1C2partialAll(temp1);

        double[] partialG2_all = temp2; // new double[4];
        for (int i = 0; i < 3; ++i) {
            double partialC1 = partialC1C2_all[i][0];
            double partialC2 = partialC1C2_all[i][1];
            double partialG2 = g1(t) * ((partialC1 * (1 + C2) + partialC2 * C1)) -
                    (partialC1 * t * expC1t * (C2 - 1) + partialC2 * (1 - expC1t)) * C1 * (1 + C2);
            double G1 = g1(t);
            partialG2 = -2 * partialG2 / (G1 * G1);
            partialG2 += partialC1;
            partialG2_all[i] = partialG2;
        }
        partialG2_all[3] = 0; // w.r.t. rho
        return partialG2_all;
    }

    // (lambda, mu, psi, rho)
    public double[] partialP0partialAll(double t, double expC1t) {
        double[] partialG2_all = partialG2partialAll(t, expC1t);

        double[] partialP0_all = temp2; // new double[4];

        double lambda = lambda();
        double G2 = g2(t);
        double mu = mu();
        double psi = psi();
//        double[] constants = getConstants();
        double C1 = getC1();
        double C2 = getC2();
        double G1 = g1(t);

//        double expC1t = Math.exp(-C1 * t); // TODO Notice this is (1) shared in many functions and (2) slow to compute

        // lambda
        partialP0_all[0] = (-mu - psi + lambda * partialG2_all[0] - G2) / (2 * lambda*lambda);
        // mu
        partialP0_all[1] = (1 + partialG2_all[1]) / (2 * lambda);
        // psi
        partialP0_all[2] = (1 + partialG2_all[2]) / (2 * lambda);
        // rho
        partialP0_all[3] = -C1 / lambda * (2 * lambda / C1 * (G1 - (1 - expC1t) * (1 + C2))) / (G1 * G1);

        return partialP0_all;
    }


    // (lambda, mu, psi, rho)
    public double[] getAllGradient(Tree tree, NodeRef node) {
        return null;
//        return getGradientLogDensityImpl((TreeModel) tree);
    }


    @Override
    public void precomputeGradientConstants() {
        precomputeConstants();
        this.savedQ = Double.MIN_VALUE;
//        this.savedPartialQ = null;
        this.partialQKnown = false;
    }

    @Override
    public void processGradientModelSegmentBreakPoint(double[] gradient, int currentModelSegment, double intervalStart, double segmentIntervalEnd) {
        return;
    }

    @Override
    public void processGradientInterval(double[] gradient, int currentModelSegment, double intervalStart, double intervalEnd, int nLineages) {
        double tOld = intervalEnd;
        double tYoung = intervalStart;
        double[] partialQ_all_old = partialQpartialAll(temp2, tOld);
        double[] partialQ_all_young;
        double Q_Old = Q(tOld);
        double Q_young;
        if (this.savedQ != Double.MIN_VALUE) {
            Q_young = this.savedQ;
        }
        else {
            Q_young = Q(tYoung);
        }
        this.savedQ = Q_Old;

        if (partialQKnown) {
            partialQ_all_young = temp3;
            System.arraycopy(partialQ, 0, partialQ_all_young, 0, 4);
        } else {
            partialQ_all_young = partialQpartialAll(temp3, tYoung);
            //System.arraycopy(partialQ_all_young, 0, savedPartialQ, 0, 4);
            partialQKnown = true;
        }
        System.arraycopy(partialQ_all_old, 0, partialQ, 0, 4);

//        if (this.savedPartialQ != null) {
//            partialQ_all_young = this.savedPartialQ;
//        }
//        else {
//            partialQ_all_young = partialQpartialAll(tYoung);
//        }
//        this.savedPartialQ = partialQ_all_old;

        for (int j = 0; j < 4; ++j) {
            gradient[j] += nLineages*(partialQ_all_young[j] / Q_young - partialQ_all_old[j] / Q_Old);
        }

    }

    @Override
    public void processGradientSampling(double[] gradient, int currentModelSegment, double intervalEnd) {
        double r = r();
        double rho = rho();
        double psi = psi();
        double t = intervalEnd;

        double timeZeroTolerance = Double.MIN_VALUE;
        boolean noSamplingAtPresent = rho() < Double.MIN_VALUE;

        if (noSamplingAtPresent || t > timeZeroTolerance) {
            double expC1t = Math.exp(-getC1() * t);
            double[] partialP0_all = partialP0partialAll(t, expC1t);
            double P0 = p0(t);
            double v = (1 - r) / ((1 - r) * P0 + r);
            for (int j = 0; j < 4; ++j) {
                gradient[j] += v * partialP0_all[j];
            }
            gradient[2] += 1 / psi;
            gradient[4] += (1 - P0) / ((1 - r)*P0 + r);
        } else {
            if (rho != 0) {
                gradient[3] += 1 / rho;
            }
        }

    }

    @Override
    public void processGradientCoalescence(double[] gradient, int currentModelSegment, double intervalEnd) {
        gradient[0] += 1 / lambda();
    }

    @Override
    public void processGradientOrigin(double[] gradient, int currentModelSegment, double totalDuration) {
        double origin = originTime.getValue(0);
        double[] partialQ_all_origin = partialQpartialAll(temp2, origin);
        double[] partialQ_all_root = partialQpartialAll(temp3, totalDuration);
        for (int i = 0; i < 4; ++i) {
            // partialLL_all[i] = 1 / (1 - p0) * partialP0_all_origin[i];
            gradient[i] += partialQ_all_root[i]/Q(totalDuration) - partialQ_all_origin[i] / Q(origin);
        }

    }

    @Override
    public void logConditioningProbability(double[] gradient) {
        return;
    }

}