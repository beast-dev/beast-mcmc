/*
 * BirthDeathSamplingModel.java
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

// Incorporated into the base model

/**
 * The reconstructed birth-death model with incomplete sampling.
 */
//public class BirthDeathSamplingModel extends BirthDeathGernhard08Model {
//
//    public static final String BIRTH_DEATH_SAMPLING_MODEL = "birthDeathSamplingModel";
//
//    private Parameter samplingProportion;
//
//    public BirthDeathSamplingModel(
//            Parameter birthDiffRateParameter,
//            Parameter relativeDeathRateParameter,
//            Parameter samplingProportion,
//            Type units) {
//        this(BIRTH_DEATH_SAMPLING_MODEL, birthDiffRateParameter,
//                relativeDeathRateParameter, samplingProportion, units);
//    }
//
//
//    BirthDeathSamplingModel(String name,
//                            Parameter birthDiffRateParameter,
//                            Parameter relativeDeathRateParameter,
//                            Parameter samplingProportion,
//                            Type units) {
//
//        super(name, birthDiffRateParameter, relativeDeathRateParameter, units);
//
//        this.samplingProportion = samplingProportion;
//        addVariable(samplingProportion);
//        samplingProportion.addBounds(new Parameter.DefaultBounds(1.0, 0.0, 1));
//    }
//
//    public final double getS() {
//        return samplingProportion.getParameterValue(0);
//    }
//
//    public final double getLambda() {
//
//        double lambda = getR() / (1 - getA());
//
//        return lambda;
//    }
//
//    public final double getMu() {
//
//        double mu = getR() * getA() / (1 - getA());
//
//        return mu;
//    }
//
//
//    public double logTreeProbability(int taxonCount) {
//
//        int n = taxonCount;
//
//        return logGamma(n + 1) + (2 * n - 1) * Math.log(getR()) + Math.log(getLambda() * getS()) * (n - 1);
//
//    }
//
//    public double logNodeProbability(Tree tree, NodeRef node) {
//        final double height = tree.getNodeHeight(node);
//        final double mrh = -getR() * height;
//        final double z = Math.log(getS() * getLambda() + (getLambda() * (1 - getS()) - getMu()) * Math.exp(mrh));
//        double l = -2 * z + mrh;
//
//        if (tree.getRoot() == node) {
//            l += mrh - z;
//        }
//        return l;
//    }
//
//    public boolean includeExternalNodesInLikelihoodCalculation() {
//        return false;
//    }
//}
