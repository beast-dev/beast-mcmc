/*
 * DiscreteTraitBranchSubstitutionParameterGradient.java
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
import dr.inference.model.BranchParameter;
import dr.inference.model.Parameter;
import dr.math.NumericalDerivative;
import dr.xml.Reportable;

/**
 * @author Marc A. Suchard
 * @author Xiang Ji
 */

//TODO: combine this with HyperParameterBranchRateGradient to remove code duplication
public abstract class HyperParameterBranchSubstitutionParameterGradient extends DiscreteTraitBranchSubstitutionParameterGradient
        implements GradientWrtParameterProvider, HessianWrtParameterProvider, Reportable, Loggable {

    protected final ArbitraryBranchRates.BranchRateTransform locationScaleTransform;
    protected final Parameter hyperParameter;

    public HyperParameterBranchSubstitutionParameterGradient(String traitName,
                                                             TreeDataLikelihood treeDataLikelihood,
                                                             BeagleDataLikelihoodDelegate likelihoodDelegate,
                                                             BranchParameter branchParameter,
                                                             Parameter hyperParameter,
                                                             boolean useHessian) {

        super(traitName, treeDataLikelihood, likelihoodDelegate,  branchParameter, useHessian);

        locationScaleTransform = branchParameter.getTransform();

        this.hyperParameter = hyperParameter;
    }

    @Override
    public double[] getGradientLogDensity() { //TODO：fix commented code chunk
//        double[] result = new double[hyperParameter.getDimension()];
//        double[] nodeGradients = super.getGradientLogDensity();
//
//        for (int i = 0; i < result.length; ++i) {
//            final int nodeNum = parameterIndexHelper.getNodeNumberFromParameterIndex(i);
//            final NodeRef node = tree.getNode(nodeNum);
//            double[] hyperChainGradient = getDifferential(tree, node);
//            for (int j = 0; j < result.length; j++) {
//                result[j] += nodeGradients[i] * hyperChainGradient[j];
//            }
//        }
//        return result;

        double[] result = new double[hyperParameter.getDimension()];
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
    public double[] getDiagonalHessianLogDensity() {
        // cannot avoid calculating full hessian in this case, use numerical method for now
        // TODO: maybe add Hessian into BEAGLE ?
        return NumericalDerivative.diagonalHessian(numeric, branchParameter.getParameterValues());
    }

    protected double getChainGradient(Tree tree, NodeRef node) {
        return 1.0;
    }


    abstract double[] getDifferential(Tree tree, NodeRef node);

}
