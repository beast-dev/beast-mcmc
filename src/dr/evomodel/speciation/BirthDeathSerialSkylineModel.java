/*
 * BirthDeathSerialSkylineModel.java
 *
 * Copyright (C) 2002-2011 Alexei Drummond and Andrew Rambaut
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
 * BEAST is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
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

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Set;

/**
 * Beginning of tree prior for birth-death + serial sampling + extant sample proportion. More Tanja magic...
 *
 * @author Alexei Drummond
 */
public class BirthDeathSerialSkylineModel extends SpeciationModel {

    // times
    Variable<Double> times;

    // birth rate
    Variable<Double> lambda;

    // death rate
    Variable<Double> mu;

    // serial sampling rate
    Variable<Double> psi;

    // extant sampling proportion
    Variable<Double> p;

    //boolean death rate is relative?
    boolean relativeDeath = false;

    // boolean stating whether sampled individuals remain infectious, or become non-infectious
    boolean sampledIndividualsRemainInfectious = false;

    double finalTimeInterval = 0.0;

    public BirthDeathSerialSkylineModel(
            Variable<Double> times,
            Variable<Double> lambda,
            Variable<Double> mu,
            Variable<Double> psi,
            Variable<Double> p,
            boolean relativeDeath,
            boolean sampledIndividualsRemainInfectious,
            double finalTimeInterval,
            Type units) {

        this("birthDeathSerialSamplingModel", times, lambda, mu, psi, p, relativeDeath, sampledIndividualsRemainInfectious, finalTimeInterval, units);
    }

    public BirthDeathSerialSkylineModel(
            String modelName,
            Variable<Double> times,
            Variable<Double> lambda,
            Variable<Double> mu,
            Variable<Double> psi,
            Variable<Double> p,
            boolean relativeDeath,
            boolean sampledIndividualsRemainInfectious,
            double finalTimeInterval,
            Type units) {

        super(modelName, units);

        int size = times.getSize();

        if (lambda.getSize() != 1 && lambda.getSize() != size)
            throw new RuntimeException("Length of Lambda parameter should be one or equal to the size of time parameter (size = " + size + ")");

        if (mu.getSize() != 1 && mu.getSize() != size)
            throw new RuntimeException("Length of mu parameter should be one or equal to the size of time parameter (size = " + size + ")");

        this.lambda = lambda;
        addVariable(lambda);
        lambda.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, 1));

        this.mu = mu;
        addVariable(mu);
        mu.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, 1));

        this.p = p;
        addVariable(p);
        p.addBounds(new Parameter.DefaultBounds(1.0, 0.0, 1));

        this.psi = psi;
        addVariable(psi);
        psi.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, 1));

        this.relativeDeath = relativeDeath;

        this.finalTimeInterval = finalTimeInterval;

        this.sampledIndividualsRemainInfectious = sampledIndividualsRemainInfectious;

    }

    public static double p0(double b, double d, double p, double psi, double t) {
        double c1 = c1(b, d, psi);
        double c2 = c2(b, d, p, psi);

        double expc1trc2 = Math.exp(-c1 * t) * (1.0 - c2);

        return (b + d + psi + c1 * ((expc1trc2 - (1.0 + c2)) / (expc1trc2 + (1.0 + c2)))) / (2.0 * b);
    }

    public static double q(double b, double d, double p, double psi, double t) {
        double c1 = c1(b, d, psi);
        double c2 = c2(b, d, p, psi);
        double res = 2 * (1 - c2 * c2) + Math.exp(-c1 * t) * (1 - c2) * (1 - c2) + Math.exp(c1 * t) * (1 + c2) * (1 + c2);
        return res;
    }

    public double p0(double t) {
        return p0(birth(), death(), p(), psi(), t);
    }

    public double q(double t) {
        return q(birth(), death(), p(), psi(), t);
    }

    private static double c1(double b, double d, double psi) {
        return Math.abs(Math.sqrt(Math.pow(b - d - psi, 2.0) + 4.0 * b * psi));
    }

    private static double c2(double b, double d, double p, double psi) {
        return -(b - d - 2.0 * b * p - psi) / c1(b, d, psi);
    }

    private double c1() {
        return c1(birth(), death(), psi());
    }

    private double c2() {
        return c2(birth(), death(), p(), psi());
    }

    public double birth() {
        return lambda.getValue(0);
    }

    public double death() {
        return relativeDeath ? mu.getValue(0) * birth() : mu.getValue(0);
    }

    public double psi() {
        return psi.getValue(0);
    }

    public double p() {

        if (finalTimeInterval == 0.0) return p.getValue(0);
        return 0;
    }

    /**
     * @param t
     * @return the birth parameter for the given time
     */
    public double lambda(double t) {

        List<Double> endTime = getEndTimes();
        int epoch = Collections.binarySearch(endTime, t);

        if (epoch < 0) {
            epoch = -epoch - 1;

        }
        return lambda.getValue(epoch);

//        int index = 0;
//        while (time > times.getValue(index)) {
//            index += 1;
//        }
//        return mu.getValue(index);
    }

    /**
     * @param t
     * @return the mutation parameter for the given time
     */
    public double mu(double t) {

        List<Double> endTime = getEndTimes();

        int epoch = Collections.binarySearch(endTime, t);

        if (epoch < 0) {
            epoch = -epoch - 1;

        }
        return mu.getValue(epoch);

//        int index = 0;
//        while (time > times.getValue(index)) {
//            index += 1;
//        }
//        return mu.getValue(index);
    }

    /**
     * Generic likelihood calculation
     *
     * @param tree the tree to calculate likelihood of
     * @return log-likelihood of density
     */
    public final double calculateTreeLogLikelihood(Tree tree) {

        //System.out.println("calculating tree log likelihood");
        double time = finalTimeInterval;

        // extant leaves
        int n = 0;
        // extinct leaves
        int m = 0;

        for (int i = 0; i < tree.getExternalNodeCount(); i++) {
            NodeRef node = tree.getExternalNode(i);
            if (tree.getNodeHeight(node) + time == 0.0) {
                n += 1;
            } else {
                m += 1;
            }
        }

        double x1 = tree.getNodeHeight(tree.getRoot()) + time;
        double c1 = c1();
        double c2 = c2();
        double b = birth();
        double p = p();

        double bottom = c1 * (c2 + 1) * (1 - c2 + (1 + c2) * Math.exp(c1 * x1));
        double logL = Math.log(1 / bottom);

        if (n > 0) {
            logL += n * Math.log(4 * p);
        }
        for (int i = 0; i < tree.getInternalNodeCount(); i++) {
            double x = tree.getNodeHeight(tree.getInternalNode(i)) + time;


            logL += Math.log(b / q(x));
        }


        for (int i = 0; i < tree.getExternalNodeCount(); i++) {
            double y = tree.getNodeHeight(tree.getExternalNode(i)) + time;

            if (y > 0.0) {

                if (sampledIndividualsRemainInfectious) { // i.e. modification (i) or (ii)
                    logL += Math.log(psi() * q(y) * p0(y));
                } else {

                    logL += Math.log(psi() * q(y));
                }
            }
        }

        return logL;
    }

    public double calculateTreeLogLikelihood(Tree tree, Set<Taxon> exclude) {
        if (exclude.size() == 0) return calculateTreeLogLikelihood(tree);
        throw new RuntimeException("Not implemented!");
    }

    public List<Double> getEndTimes() {
        List<Double> endTimes = new ArrayList<Double>();
        for (int i = 0; i < endTimes.size(); i++) {
            endTimes.add(times.getValue(i));
        }
        return endTimes;
    }

    public static void main(String[] args) {

        // test framework

        Variable<Double> times = new Variable.D(1, 10);

        Variable<Double> mu = new Variable.D(1, 10);
        for (int i = 0; i < mu.getSize(); i++) {
            mu.setValue(i, i + 1.0);
        }

        Variable<Double> lambda = new Variable.D(1, 10);
        Variable<Double> psi = new Variable.D(0.5, 1);
        Variable<Double> p = new Variable.D(0.5, 1);
        boolean relativeDeath = false;
        boolean sampledIndividualsRemainInfectious = false;
        double finalTime = 50;

        BirthDeathSerialSkylineModel model =
                new BirthDeathSerialSkylineModel(times, lambda, mu, psi, p, relativeDeath,
                        sampledIndividualsRemainInfectious, finalTime, Type.SUBSTITUTIONS);

        for (double time = 0.5; time < 10; time += 1) {
            System.out.println("mu at time " + time + " is " + model.mu(time));
        }

    }
}