/*
 * BirthDeathGernhard08Model.java
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
import dr.inference.model.Parameter;
import dr.inference.model.Variable;

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

        return birth() + death() + psi() + c1 * ((expc1trc2 - 1.0) / (expc1trc2 + 1.0));
    }

    private double c1() {

        double b = birth();

        return Math.abs(Math.sqrt(Math.pow(b - death() - psi(), 2) + 4 * b * psi()));
    }

    private double c2() {

        double b = birth();

        return -(b - death() - 2 * b * p() - psi()) / c1();
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

    public double logTreeProbability(int taxonCount) {
        throw new RuntimeException("to be implemented");
    }

    public double logNodeProbability(Tree tree, NodeRef node) {
        throw new RuntimeException("to be implemented");
    }

    public boolean includeExternalNodesInLikelihoodCalculation() {
        return false;
    }
}