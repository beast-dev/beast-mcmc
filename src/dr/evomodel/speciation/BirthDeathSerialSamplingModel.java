/*
 * BirthDeathSerialSamplingModel.java
 *
 * Copyright (C) 2002-2009 Alexei Drummond and Andrew Rambaut
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

import java.util.Set;

/**
 * Beginning of tree prior for birth-death + serial sampling + extant sample proportion. More Tanja magic...
 *
 * @author Alexei Drummond
 */
public class BirthDeathSerialSamplingModel extends SpeciationModel {

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

    public BirthDeathSerialSamplingModel(
            Variable<Double> lambda,
            Variable<Double> mu,
            Variable<Double> psi,
            Variable<Double> p,
            boolean relativeDeath,
            Type units) {

        this("birthDeathSerialSamplingModel", lambda, mu, psi, p, relativeDeath, units);
    }

    public BirthDeathSerialSamplingModel(
            String modelName,
            Variable<Double> lambda,
            Variable<Double> mu,
            Variable<Double> psi,
            Variable<Double> p,
            boolean relativeDeath,
            Type units) {

        super(modelName, units);

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

    }


    public static double p0(double b, double d, double p, double psi, double t) {
        double c1 = c1(b, d, psi);
        double c2 = c2(b, d, p, psi);

        double expc1trc2 = Math.exp(-c1 * t) * ((1.0 - c2) / (1.0 + c2));

        return (b + d + psi + c1 * ((expc1trc2 - 1.0) / (expc1trc2 + 1.0))) / (2.0 * b);
    }

    public static double p1(double b, double d, double p, double psi, double t) {
        double c1 = c1(b, d, psi);
        double c2 = c2(b, d, p, psi);
        double c3 = c3(b, d, p, psi);
        double pc2 = 1.0 + c2;
        double mc2 = 1.0 - c2;
        double numerator = -4.0 * p * c1 * c1;
        double denominator = c3 * (2.0 + Math.exp(-c1 * t) * (mc2 / pc2) + Math.exp(c1 * t) * (pc2 / mc2));
        return numerator / denominator;
    }

    public double p0(double t) {
        return p0(birth(), death(), p(), psi(), t);
    }

    public double p1(double t) {
        return p1(birth(), death(), p(), psi(), t);
    }

    public double q(double s, double t) {
        return psi() * p0(s) * p1(t) / p1(s);
    }

    private static double c1(double b, double d, double ss) {
        return Math.abs(Math.sqrt(Math.pow(b - d - ss, 2.0) + 4.0 * b * ss));
    }

    private static double c2(double b, double d, double p, double ss) {
        return -(b - d - 2.0 * b * p - ss) / c1(b, d, ss);
    }

    private static double c3(double b, double d, double p, double ss) {
        return 4.0 * b * (p * (d + b * (p - 1.0) + ss) - ss);
    }


    private double c1() {
        return c1(birth(), death(), psi());
    }

    private double c2() {
        return c2(birth(), death(), p(), psi());
    }

    private double c3() {
        return c3(birth(), death(), p(), psi());
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
        return p.getValue(0);
    }

    /**
     * Generic likelihood calculation
     *
     * @param tree the tree to calculate likelihood of
     * @return log-likelihood of density
     */
    public final double calculateTreeLogLikelihood(Tree tree) {

        //System.out.println("calculating tree log likelihood");

        // extant leaves
        int n = 0;

        // extinct leaves
        int m = 0;

        for (int i = 0; i < tree.getExternalNodeCount(); i++) {
            NodeRef node = tree.getExternalNode(i);
            if (tree.getNodeHeight(node) == 0.0) {
                n += 1;
            } else {
                m += 1;
            }
        }

        //System.out.println("m = " + m);
        //System.out.println("n = " + n);

        double x1 = tree.getNodeHeight(tree.getRoot());
        double c1 = c1();
        double c2 = c2();
        double c3 = c3();
        double b = birth();


        double top = 4.0 * c1 * (c2 - 1.0) * p();
        double bottom = c3 * (1.0 - c2 + (1.0 + c2) * Math.exp(c1 * x1));
        //System.out.println("top=" + top);
        //System.out.println("bottom=" + bottom);

        double logL = Math.log(n * b * top / bottom);
        //System.out.println("logL=" + logL);
        for (int i = 0; i < tree.getInternalNodeCount(); i++) {
            double x = tree.getNodeHeight(tree.getInternalNode(i));

            logL += Math.log(b * p1(x));
        }
        for (int i = 0; i < tree.getExternalNodeCount(); i++) {
            double y = tree.getNodeHeight(tree.getExternalNode(i));

            if (y > 0.0) {
                logL += Math.log(psi() * p0(y) / p1(y));
            }
        }
        return logL;
    }

    public double calculateTreeLogLikelihood(Tree tree, Set<Taxon> exclude) {
        if (exclude.size() == 0) return calculateTreeLogLikelihood(tree);
        throw new RuntimeException("Not implemented!");
    }
}