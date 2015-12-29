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
import dr.evolution.tree.Tree;
import dr.evolution.util.Taxon;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;

import java.io.IOException;
import java.util.Arrays;
import java.util.Set;

/**
 * Beginning of tree prior for birth-death + serial sampling + extant sample proportion. More Tanja magic...
 * <p/>
 * log:
 * 25 Mar 2011, Denise: added int i (index) for the Variables that change over time such as the methods p0(..), q(..); fixed some formulas (old versions commented out)
 * unclear marked with todo
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

    Variable<Double> origin;

    //boolean death rate is relative?
    boolean relativeDeath = false;

    // the number of intervals;
    int size = 1;

    double t_root;
    double x0;
    protected double[] p0_iMinus1;
    protected double[] Ai;
    protected double[] Bi;

    protected boolean birthChanges = true;
    protected boolean deathChanges = true;
    protected boolean samplingChanges = true;

    protected boolean timesStartFromOrigin = true;
    protected double[] timesFromTips;

    public BirthDeathSerialSkylineModel(
            Variable<Double> times,
            Variable<Double> lambda,
            Variable<Double> mu,
            Variable<Double> psi,
            Variable<Double> p,
            Variable<Double> origin,
            boolean relativeDeath,
            boolean sampledIndividualsRemainInfectious,
            boolean timesStartFromOrigin,
            Type units) {

        this("birthDeathSerialSamplingModel", times, lambda, mu, psi, p, origin, relativeDeath,
                sampledIndividualsRemainInfectious, timesStartFromOrigin, units);
    }

    public BirthDeathSerialSkylineModel(
            String modelName,
            Variable<Double> times,
            Variable<Double> lambda,
            Variable<Double> mu,
            Variable<Double> psi,
            Variable<Double> p,
            Variable<Double> origin,
            boolean relativeDeath,
            boolean sampledIndividualsRemainInfectious,
            boolean timesStartFromOrigin,
            Type units) {

        super(modelName, units);

        this.size = times.getSize();

        if (lambda.getSize() != 1 && lambda.getSize() != size)
            throw new RuntimeException("Length of Lambda parameter should be one or equal to the size of time parameter (size = " + size + ")");

        if (mu.getSize() != 1 && mu.getSize() != size)
            throw new RuntimeException("Length of mu parameter should be one or equal to the size of time parameter (size = " + size + ")");

        this.timesStartFromOrigin = timesStartFromOrigin;

        this.times = times;
        addVariable(times);
        times.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, times.getSize()));

        this.lambda = lambda;
        addVariable(lambda);
        lambda.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, lambda.getSize()));

        this.mu = mu;
        addVariable(mu);
        mu.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, mu.getSize()));

        this.p = p;
        addVariable(p);
        p.addBounds(new Parameter.DefaultBounds(1.0, 0.0, p.getSize()));

        this.origin = origin;
        addVariable(origin);
        p.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, origin.getSize()));

        this.psi = psi;
        addVariable(psi);
        psi.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, psi.getSize()));


        this.relativeDeath = relativeDeath;
    }

    /**
     * @param time the time
     * @param tree the tree
     * @return the number of lineages that exist at the given time in the given tree.
     */
    public int lineageCountAtTime(double time, Tree tree) {

        int count = 1;
        for (int i = 0; i < tree.getInternalNodeCount(); i++) {
            if (tree.getNodeHeight(tree.getInternalNode(i)) > time) count += 1;
        }
        for (int i = 0; i < tree.getExternalNodeCount(); i++) {
            if (tree.getNodeHeight(tree.getExternalNode(i)) > time) count -= 1;
        }
        return count;
    }

    public double Ai(double b, double g, double psi) {

        return Math.sqrt((b - g - psi) * (b - g - psi) + 4.0 * b * psi);
    }

    public double Bi(double b, double g, double psi, double A, double p0) {

        return (-((1.0 - 2.0 * p0) * b + g + psi) / A);
    }

    public double p0(int index, double t, double ti) {

        return p0(birth(birthChanges ? index : 0), death(deathChanges ? index : 0), psi(samplingChanges ? index : 0), Ai[index], Bi[index], t, ti);
    }

    public double p0(double b, double g, double psi, double A, double B, double t, double ti) {

        return ((b + g + psi - A * ((Math.exp(A * (t - ti)) * (1.0 - B) - (1.0 + B))) / (Math.exp(A * (t - ti)) * (1.0 - B) + (1.0 + B))) / (2.0 * b));
    }

    public double g(int index, double t, double ti) {

        double oneMinusBiSq = (1.0 - Bi[index]) * (1.0 - Bi[index]);
        double onePlusBiSq = (1.0 + Bi[index]) * (1.0 + Bi[index]);

        return 4.0 / (2.0 * (1.0 - Bi[index] * Bi[index]) + Math.exp(Ai[index] * (t - ti)) *
                oneMinusBiSq + Math.exp(-Ai[index] * (t - ti)) * onePlusBiSq);
    }

    /**
     * Returns the time at which epoch i begins. If
     *
     * @param i index of the epoch
     * @return the time at which this epoch begins
     */
    public double t(int i) {
        return timesFromTips[i];
    }

    public double birth(int i) {
        return lambda.getValue(i);
    }

    public double death(int i) {
        return relativeDeath ? mu.getValue(i) * birth(i) : mu.getValue(i);
    }

    public double psi(int i) {
        return psi.getValue(i);
    }

    public double p() {
        return p.getValue(0);
    }

    /**
     * @param t
     * @return the birth parameter for the given time
     */
    public double lambda(double t) {

        return lambda.getValue(index(t));

    }

    /**
     * @param t
     * @return the mutation parameter for the given time
     */
    public double mu(double t) {

        return mu.getValue(index(t));
    }

    public int index(double t) {

        int epoch = Arrays.binarySearch(timesFromTips, t);

        if (epoch < 0) {
            epoch = -epoch - 1;
        }
        return Math.max(epoch - 1, 0);
    }

    /*    calculate and store Ai, Bi and p0_iMinus1        */
    public void preCalculation(Tree tree) {

        t_root = tree.getNodeHeight(tree.getRoot());
        x0 = t_root + origin.getValue(0);

        // set up timesFromTips array
        if (timesFromTips == null) {
            timesFromTips = new double[times.getSize()];
        }
        if (timesStartFromOrigin) {
            timesFromTips[0] = 0;
            for (int i = 1; i < timesFromTips.length; i++) {
                timesFromTips[i] = Math.max(0, x0 - times.getValue(timesFromTips.length - i));
            }

        } else {
            for (int i = 0; i < timesFromTips.length; i++) {
                timesFromTips[i] = times.getValue(i);
            }
        }

        Ai = new double[size];
        Bi = new double[size];
        p0_iMinus1 = new double[size];

        for (int i = 0; i < size; i++) {
            Ai[i] = Ai(birth(birthChanges ? i : 0), death(deathChanges ? i : 0), psi(samplingChanges ? i : 0));

            //System.out.println("Ai[" + i + "]=" + Ai[i]);
        }

        Bi[0] = Bi(birth(0), death(0), psi(0), Ai[0], 1);
        //System.out.println("Bi[0]=" + Bi[0]);
        for (int i = 1; i < size; i++) {
            p0_iMinus1[i - 1] = p0(birth(birthChanges ? (i - 1) : 0), death(deathChanges ? (i - 1) : 0), psi(samplingChanges ? (i - 1) : 0), Ai[i - 1], Bi[i - 1], t(i), t(i - 1));
            Bi[i] = Bi(birth(birthChanges ? i : 0), death(deathChanges ? i : 0), psi(samplingChanges ? i : 0), Ai[i], p0_iMinus1[i - 1]);

            //System.out.println("Bi[" + i + "]=" + Bi[i]);
        }
    }

    /**
     * Generic likelihood calculation
     *
     * @param tree the tree to calculate likelihood of
     * @return log-likelihood of density
     */
    public final double calculateTreeLogLikelihood(Tree tree) {

        // number of lineages at each time ti
        int[] n = new int[size];
        int nTips = tree.getExternalNodeCount();
        preCalculation(tree);

        int index = size - 1;      // x0 must be in last interval

        double t = t(index);
        double g = g(index, x0, t);

        double logP = Math.log(g);

        // first product term in f[T]
        for (int i = 0; i < tree.getInternalNodeCount(); i++) {

            double x = tree.getNodeHeight(tree.getInternalNode(i));
            index = index(x);

            double contrib = Math.log(birth(birthChanges ? index : 0) * g(index, x, t(index)));

            logP += contrib;

            //System.out.println("internalNode.logP=" + contrib);

            t = t(index);
            g = g(index, x, t);


            //System.out.println("logP+=" + (Math.log(birth(birthChanges ? index : 0) * g)) + " t= " + t + " g=" + g);
        }

        // middle product term in f[T]
        for (int i = 0; i < nTips; i++) {

            double y = tree.getNodeHeight(tree.getExternalNode(i));
            index = index(y);

            double contrib = Math.log(psi(samplingChanges ? index : 0)) - Math.log(g(index, y, t(index)));
            ;
            logP += contrib;

            //System.out.println("externalNode.logP=" + contrib);

        }

        // last product term in f[T], factorizing from 1 to m
        for (int j = 0; j < size - 1; j++) {

            double contrib = 0;

            double time = t(j + 1);
            n[j] = lineageCountAtTime(time, tree);
            if (n[j] > 0) {
                contrib += n[j] * Math.log(g(j, time, t(j)));
                //System.out.println("n[" + j + "]" + n[j] + " time=" + time + " t(" + j + ")=" + t(j));
            }

            logP += contrib;
            //System.out.println("last term=" + contrib);
        }
        return logP;
    }

    public double calculateTreeLogLikelihood(Tree tree, Set<Taxon> exclude) {
        if (exclude.size() == 0) return calculateTreeLogLikelihood(tree);
        throw new RuntimeException("Not implemented!");
    }

    public static void main(String[] args) throws IOException, Importer.ImportException {

        // test framework

        Variable<Double> times = new Variable.D(1, 10);

        Variable<Double> mu = new Variable.D(1, 10);
        for (int i = 0; i < mu.getSize(); i++) {
            times.setValue(i, (i + 1) * 2.0);
            mu.setValue(i, i + 1.0);
        }

        Variable<Double> lambda = new Variable.D(1, 10);
        Variable<Double> psi = new Variable.D(0.5, 1);
        Variable<Double> p = new Variable.D(0.5, 1);
        Variable<Double> origin = new Variable.D(0.5, 1);
        boolean relativeDeath = false;
        boolean sampledIndividualsRemainInfectious = false;
        boolean timesStartFromOrigin = false;

        BirthDeathSerialSkylineModel model =
                new BirthDeathSerialSkylineModel(times, lambda, mu, psi, p, origin, relativeDeath,
                        sampledIndividualsRemainInfectious, timesStartFromOrigin, Type.SUBSTITUTIONS);

        NewickImporter importer = new NewickImporter("((A:6,B:5):4,(C:3,D:2):1);");
        Tree tree = importer.importNextTree();

        model.calculateTreeLogLikelihood(tree);

        for (int i = 0; i < times.getSize(); i += 1) {
            System.out.println("mu at time " + i + " is " + model.mu(i));
            System.out.println("p0 at time " + i + " is " + model.p0(0, i, i));
        }

    }
}