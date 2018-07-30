/*
 * LocationScaleGradient.java
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

package dr.evomodel.treedatalikelihood.discrete;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evomodel.branchratemodel.ArbitraryBranchRates;
import dr.evomodel.treedatalikelihood.BeagleDataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.inference.hmc.GradientWrtParameterProvider;
import dr.inference.hmc.HessianWrtParameterProvider;
import dr.inference.loggers.Loggable;
import dr.inference.model.CompoundParameter;
import dr.inference.model.Parameter;
import dr.xml.Reportable;

/**
 * @author Marc A. Suchard
 * @author Xiang Ji
 */

public abstract class HyperParameterBranchRateGradient extends DiscreteTraitBranchRateGradient
        implements GradientWrtParameterProvider, HessianWrtParameterProvider, Reportable, Loggable {

    protected final ArbitraryBranchRates.BranchRateTransform.LocationScaleLogNormal locationScaleTransform;

    protected HyperParameterBranchRateGradient(String traitName,
                                             TreeDataLikelihood treeDataLikelihood,
                                             BeagleDataLikelihoodDelegate likelihoodDelegate,
                                             Parameter parameter,
                                             boolean useHessian) {

        super(traitName, treeDataLikelihood, likelihoodDelegate, parameter, useHessian);

        if (!(branchRateModel.getTransform() instanceof ArbitraryBranchRates.BranchRateTransform.LocationScaleLogNormal)) {
            throw new IllegalArgumentException("Must provide a LocationScaleLogNormal transform.");
        }

        locationScaleTransform = (ArbitraryBranchRates.BranchRateTransform.LocationScaleLogNormal) branchRateModel.getTransform();
    }

    @Override
    public double[] getGradientLogDensity() {
        double[] result = new double[rateParameter.getDimension()];
        double[] nodeGradients = super.getGradientLogDensity();
        int v = 0;
        for (int i = 0; i < tree.getNodeCount(); ++i) {
            final NodeRef node = tree.getNode(i);
            if (!tree.isRoot(node)) {
                double[] hyperChainGradient = getDifferential(tree, node);
                for (int j = 0; j < result.length; j++) {
                    result[j] += nodeGradients[v] * hyperChainGradient[j];
                }
                v++;
            }
        }
        return result;
    }

    @Override
//    public double[] getDiagonalHessianLogDensity() {  // TODO: fix it
//        double[] result = new double[rateParameter.getDimension()];
//        double[] nodeGradients = super.getGradientLogDensity();
//        double[] nodeDiagonalHessian = super.getDiagonalHessianLogDensity();
//        int v = 0;
//        for (int i = 0; i < tree.getNodeCount(); ++i) {
//            final NodeRef node = tree.getNode(i);
//            if (!tree.isRoot(node)) {
//                double[] hyperChainGradient = getDifferential(tree, node);
//                double[] hyperChainDiagonalHessian = getSecondDifferential(tree, node);
//                for (int j = 0; j < result.length; j++) {
//                    result[j] += nodeGradients[v] * hyperChainDiagonalHessian[j] + nodeDiagonalHessian[v] * hyperChainGradient[j] * hyperChainGradient[j];
////                    result[j] += nodeGradients[v] * hyperChainGradient[j];
//                }
//                v++;
//            }
//        }
//        return result;
//    }

    public double[] getDiagonalHessianLogDensity() {

        double[] result = new double[rateParameter.getDimension()];

        //Do single call to traitProvider with node == null (get full tree)
        double[] diagonalHessian = (double[]) treeDataLikelihood.getTreeTrait("Hessian").getTrait(tree, null);
        double[] gradient = (double[]) treeTraitProvider.getTrait(tree, null);

        int v = 0;
        for (int i = 0; i < tree.getNodeCount(); ++i) {
            final NodeRef node = tree.getNode(i);
            if (!tree.isRoot(node)) {
                final double chainGradient = getChainGradient(tree, node);
                double[] hyperChainGradient = getDifferential(tree, node);
                double[] hyperChainDiagonalHessian = getSecondDifferential(tree, node);
                for (int j = 0; j < result.length; j++) {
                    final double jointGradient = chainGradient * hyperChainGradient[j];
                    result[j] += diagonalHessian[v] * jointGradient * jointGradient + gradient[v] * chainGradient * hyperChainDiagonalHessian[j];
                }
                v++;
            }
        }

        return result;
    }


    abstract double[] getDifferential(Tree tree, NodeRef node);

    abstract double[] getSecondDifferential(Tree tree, NodeRef node);

}
