/*
 * NewHamiltonianMonteCarloOperator.java
 *
 * Copyright (c) 2002-2017 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.inference.operators.hmc;

import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeTrait;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.evomodel.treedatalikelihood.continuous.ContinuousDataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.continuous.MultivariateTraitDebugUtilities;
import dr.evomodel.treedatalikelihood.preorder.TipFullConditionalDistributionDelegate;
import dr.evomodel.treedatalikelihood.preorder.TipGradientViaFullConditionalDelegate;
import dr.inference.hmc.GradientWrtParameterProvider;
import dr.inference.model.MatrixParameterInterface;
import dr.inference.model.Parameter;
import dr.inference.operators.AbstractCoercableOperator;
import dr.inference.operators.CoercionMode;
import dr.inference.operators.SimpleMCMCOperator;
import dr.math.distributions.NormalDistribution;
import dr.util.Transform;

/**
 * @author Zhenyu Zhang
 * @author Marc A. Suchard
 */

public class NewBouncyParticleOperator extends SimpleMCMCOperator {

//    final GradientWrtParameterProvider gradientProvider;
    final MatrixParameterInterface precisionMatrix;
    // Something that returns the conditional distribution
    final double[][] sigma0;

    private final TreeDataLikelihood treeDataLikelihood;
    private final ContinuousDataLikelihoodDelegate likelihoodDelegate;

    private final TreeTrait gradientProvider;
    private final TreeTrait densityProvider;
    private final Tree tree;
    
    public NewBouncyParticleOperator(CoercionMode mode, double weight,
                                     TreeDataLikelihood treeDataLikelihood,
                                     ContinuousDataLikelihoodDelegate likelihoodDelegate,
                                     String traitName) {

        setWeight(weight);

        this.treeDataLikelihood  = treeDataLikelihood;
        this.likelihoodDelegate = likelihoodDelegate;

        String gradientName = TipGradientViaFullConditionalDelegate.getName(traitName);
        if (treeDataLikelihood.getTreeTrait(gradientName) == null) {
            likelihoodDelegate.addFullConditionalGradientTrait(traitName);
        }
        gradientProvider = treeDataLikelihood.getTreeTrait(gradientName);

        assert (gradientProvider != null);

        String fcdName = TipFullConditionalDistributionDelegate.getName(traitName);
        if (treeDataLikelihood.getTreeTrait(fcdName) == null) {
            likelihoodDelegate.addFullConditionalDensityTrait(traitName);
        }
        densityProvider = treeDataLikelihood.getTreeTrait(fcdName);

        assert (densityProvider != null);

        this.tree = treeDataLikelihood.getTree();
        this.precisionMatrix = likelihoodDelegate.getPrecisionParameter();

        this.sigma0 = getTreeVariance();

    }

    private double[][] getTreeVariance() {

        double priorSampleSize = likelihoodDelegate.getRootProcessDelegate().getPseudoObservations();

        return MultivariateTraitDebugUtilities.getTreeVariance(tree, 1.0,
                /*Double.POSITIVE_INFINITY*/ priorSampleSize);
    }

    @Override
    public String getPerformanceSuggestion() {
        return null;
    }

    @Override
    public String getOperatorName() {
        return "Bouncy Particle operator";
    }

    @Override
    public double doOperation() {


        // Get all gradients

//        for (int taxon = 0; taxon < nTaxa; ++taxon) {
//            double[] taxonGradient = (double[]) gradientProvider.getTrait(tree, tree.getExternalNode(taxon));
//            System.arraycopy(taxonGradient, 0, gradient, offsetOutput, taxonGradient.length);
//            offsetOutput += taxonGradient.length;
//        }

//        return leapFrog();
        return 0.0;
    }

    private long count = 0;

    private static final boolean DEBUG = false;

//    private double leapFrog() {

//        if (DEBUG) {
//            if (count % 5 == 0) {
//                System.err.println(stepSize);
//            }
//            ++count;
//        }
//
//        final int dim = gradientProvider.getDimension();
//
//        final double sigmaSquared = drawDistribution.getSD() * drawDistribution.getSD();
//
//        final double[] momentum = drawInitialMomentum(drawDistribution, dim);
//        final double[] position = leapFrogEngine.getInitialPosition();
//
//        final double prop = getScaledDotProduct(momentum, sigmaSquared) +
//                leapFrogEngine.getParameterLogJacobian();
//
//        leapFrogEngine.updateMomentum(position, momentum,
//                gradientProvider.getGradientLogDensity(), stepSize / 2);
//
//
//        if (DEBUG) {
//            System.err.println("nSteps = " + nSteps);
//        }
//
//        for (int i = 0; i < nSteps; i++) { // Leap-frog
//
//            leapFrogEngine.updatePosition(position, momentum, stepSize, sigmaSquared);
//
//            if (i < (nSteps - 1)) {
//                leapFrogEngine.updateMomentum(position, momentum,
//                        gradientProvider.getGradientLogDensity(), stepSize);
//            }
//        }
//
//        leapFrogEngine.updateMomentum(position, momentum,
//                gradientProvider.getGradientLogDensity(), stepSize / 2);
//
//        final double res = getScaledDotProduct(momentum, sigmaSquared) +
//                leapFrogEngine.getParameterLogJacobian();
//
//        return prop - res; //hasting ratio
//    }

//    @Override
//    public double getCoercableParameter() {
//        return Math.log(stepSize);
//    }
//
//    @Override
//    public void setCoercableParameter(double value) {
//        if (DEBUG) {
//            System.err.println("Setting coercable paramter: " + getCoercableParameter() + " -> " + value);
//        }
//        stepSize = Math.exp(value);
//    }
//
//    @Override
//    public double getRawParameter() {
//        return 0;
//    }
}
