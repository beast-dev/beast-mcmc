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
import dr.xml.Reportable;

import java.util.Collections;
import java.util.List;
import java.util.Set;

/**
 * A phylogenetic birth-death-sampling model which includes serial sampling, sampling at present, and the possibility of treatmentProbability.
 */
public class NewBirthDeathSerialSamplingModel extends MaskableSpeciationModel implements Citable {

    // extant sampling proportion
    Variable<Double> samplingFractionAtPresent;

    // birth rate
    Variable<Double> birthRate;

    // death rate
    Variable<Double> deathRate;

    // serial sampling rate
    Variable<Double> serialSamplingRate;

    // "treatmentProbability" parameter aka r aka Pr(death | lineage is sampled)
    Variable<Double> treatmentProbability;

    // the originTime of the infection, origin > tree.getRoot();
    Variable<Double> originTime;

    private boolean conditionOnSurvival;

    // useful constants we don't want to compute nTaxa times
    private double storedC1 = Double.NEGATIVE_INFINITY;
    private double storedC2 = Double.NEGATIVE_INFINITY;

    public NewBirthDeathSerialSamplingModel(
            Variable<Double> birthRate,
            Variable<Double> deathRate,
            Variable<Double> serialSamplingRate,
            Variable<Double> treatmentProbability,
            Variable<Double> samplingFractionAtPresent,
            Variable<Double> originTime,
            boolean condition,
            Type units) {

        this("NewBirthDeathSerialSamplingModel", birthRate, deathRate, serialSamplingRate, treatmentProbability, samplingFractionAtPresent, originTime, condition, units);
    }

    public NewBirthDeathSerialSamplingModel(
            String modelName,
            Variable<Double> birthRate,
            Variable<Double> deathRate,
            Variable<Double> serialSamplingRate,
            Variable<Double> treatmentProbability,
            Variable<Double> samplingFractionAtPresent,
            Variable<Double> originTime,
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
        // TODO off by a factor of 4 from Gavryushkina et al (2014) and Magee and Hoehna (2021)
        // Should only be a problem if tree is being inferred and sampled ancestors are allowed
        double res = c1 * t + 2.0 * Math.log( Math.exp(-c1 * t) * (1.0 - c2) + (1.0 + c2) ); // operate directly in logspace, c1 * t too big
        return res;
    }

    private static double c1(double lambda, double mu, double psi) {
        return Math.abs(Math.sqrt(Math.pow(lambda - mu - psi, 2.0) + 4.0 * lambda * psi));
    }

    private static double c2(double lambda, double mu, double psi, double rho) {
        return -(lambda - 2.0 * rho * lambda - mu - psi)/c1(lambda, mu, psi);
    }

    private void precomputeConstants() {
        this.storedC1 = c1(lambda(), mu(), psi());
        this.storedC2 = c2(lambda(), mu(), psi(), rho());
    }

    public double p0(double t) {
        return p0(lambda(), mu(), psi(), rho(), storedC1, storedC2, t);
    }

    public double logq(double t) {
        return logq(storedC1, storedC2, t);
    }

    public double lambda() {
        if (mask != null) return mask.lambda();
        else {
            return birthRate.getValue(0);
        }
    }

    public double mu() {
        if (mask != null) return mask.mu();
        else {
            return deathRate.getValue(0);
        }
    }

    public double psi() {
        if (mask != null) return mask.psi();
        else {
            return serialSamplingRate.getValue(0);
        }
    }

    public double r() {
        if (mask != null) return mask.r();
        return treatmentProbability.getValue(0);
    }

    public double rho() {
        if (mask != null) return mask.rho();
        return samplingFractionAtPresent.getValue(0);
    }


    /**
     * Generic likelihood calculation
     *
     * @param tree the tree to calculate likelihood of
     * @return log-likelihood of density
     */
    public final double calculateTreeLogLikelihood(Tree tree) {
        double lambda = lambda();
        double mu = mu();
        double psi = psi();
        double r = r();
        double rho = rho();

        double timeZeroTolerance = Double.MIN_VALUE;
        boolean noSamplingAtPresent = rho < Double.MIN_VALUE;

        precomputeConstants();

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
            if (tree.getNodeHeight(node) < timeZeroTolerance) {
                n += 1;
            } else {
                m += 1;
            }
        }

        if (noSamplingAtPresent && n < 1) {
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

        // TODO tip times on starting tree are all 0?
        for (int i = 0; i < tree.getExternalNodeCount(); i++) {
            double y = tree.getNodeHeight(tree.getExternalNode(i));
            if (noSamplingAtPresent || y > timeZeroTolerance) {
                logL += Math.log(psi * (r + (1.0 - r) * p0(y))) + logq(y);
//                System.err.println("logq(y) = " + logq(y));
            }
        }
//        System.err.println("r is " + r);
//        System.err.println("new logL is " + logL);
        // TODO conditioning
        if ( conditionOnSurvival ) {
            logL -= Math.log(1.0 - p0(origin));
        }

        return logL;
    }

    public double calculateTreeLogLikelihood(Tree tree, Set<Taxon> exclude) {
        if (exclude.size() == 0) return calculateTreeLogLikelihood(tree);
        throw new RuntimeException("Not implemented!");
    }

    @Override
    public double getNodeGradient(Tree tree, NodeRef node) {
        throw new RuntimeException("Not yet implemented!");
    }

    public void mask(SpeciationModel mask) {
        if (mask instanceof NewBirthDeathSerialSamplingModel) {
            this.mask = (NewBirthDeathSerialSamplingModel) mask;
        } else {
            throw new IllegalArgumentException();
        }
    }

    public void unmask() {
        mask = null;
    }

    // if a mask exists then use the mask's parameters instead (except for originTime and finalTimeInterval)
    NewBirthDeathSerialSamplingModel mask = null;

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
}