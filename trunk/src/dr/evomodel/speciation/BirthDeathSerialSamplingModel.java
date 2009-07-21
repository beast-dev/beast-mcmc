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
 * Beginning of tree prior for birth-death + serial sampling + modern sample. More Tanja magic...
 *
 * @author Alexei Drummond
 */
public class BirthDeathSerialSamplingModel extends SpeciationModel {

    public static final String BIRTH_DEATH__SERIAL_MODEL = "birthDeathSerialSampling";

    // birth rate
    Variable<Double> lambda;

    // death rate
    Variable<Double> mu;

    // serial sampling rate
    Variable<Double> psi;

    // sampling proportion
    Variable<Double> p;

    public BirthDeathSerialSamplingModel(
            Variable<Double> lambda,
            Variable<Double> mu,
            Variable<Double> psi,
            Variable<Double> p,
            Type units) {

        this(BIRTH_DEATH__SERIAL_MODEL, lambda, mu, psi, p, units);
    }

    BirthDeathSerialSamplingModel(
            String modelName,
            Variable<Double> lambda,
            Variable<Double> mu,
            Variable<Double> psi,
            Variable<Double> p,
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

    }

    public double p0(double t) {

        double c1 = c1();
        double c2 = c2();

        double expc1trc2 = Math.exp(-c1 * t) * ((1.0 - c2) / (1.0 + c2));

        double p0 = birth() + death() + psi() + c1 * ((expc1trc2 - 1.0) / (expc1trc2 + 1.0));

        System.out.println("p0(" + t + ")=" + p0);

        return p0;
    }

    public double p1(double t) {
        double c1 = c1();
        double c2 = c2();
        double c3 = c3();
        double pc2 = 1.0 + c2;
        double mc2 = 1.0 - c2;
        double numerator = -4.0 * p() * c1 * c1;
        double denominator = c3 * (2.0 + Math.exp(-c1 * t) * (mc2 / pc2) + Math.exp(c1 * t) * (pc2 / mc2));
        return numerator / denominator;
    }

    public double q(double s, double t) {
        return psi() * p0(s) * p1(t) / p1(s);
    }

    private double c1() {

        double b = birth();
        double s = psi();

        return Math.abs(Math.sqrt(Math.pow(b - death() - s, 2.0) + 4.0 * b * s));
    }

    private double c2() {

        double b = birth();

        return -(b - death() - 2.0 * b * p() - psi()) / c1();
    }

    private double c3() {

        double b = birth();
        double s = psi();

        return 4.0 * b * (p() * (death() + b * (p() - 1.0) + s) - s);
    }


    public double birth() {
        return lambda.getValue(0);
    }

    public double death() {
        return mu.getValue(0);
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
     * @param tree
     * @return log-likelihood of density
     */
    public final double calculateTreeLogLikelihood(Tree tree) {

        System.out.println("calculating tree log likelihood");

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

        double x1 = tree.getNodeHeight(tree.getRoot());
        double c1 = c1();
        System.out.println("c1=" + c1);
        double c2 = c2();
        System.out.println("c2=" + c2);
        double c3 = c3();
        System.out.println("c3=" + c3);
        double b = birth();


        double top = 4.0 * c1 + (c2 - 1.0) * p();
        double bottom = c3 * (1.0 - c2 + (1.0 + c2) * Math.exp(c1 * x1));
        System.out.println("top=" + top);
        System.out.println("bottom=" + bottom);

        double logL = Math.log(n * b * top / bottom);
        System.out.println("logL=" + logL);
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