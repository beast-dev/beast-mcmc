/*
 * TreeTipGradient.java
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

package dr.evomodel.treedatalikelihood.continuous;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeTrait;
import dr.evomodel.branchratemodel.ArbitraryBranchRates;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.evomodel.treedatalikelihood.preorder.TipFullConditionalDistributionDelegate;
import dr.inference.hmc.GradientWrtParameterProvider;
import dr.inference.model.Likelihood;
import dr.inference.model.MatrixParameterInterface;
import dr.inference.model.Parameter;
import dr.math.matrixAlgebra.WrappedVector;
import dr.math.matrixAlgebra.missingData.MissingOps;
import dr.xml.Reportable;
import org.ejml.data.DenseMatrix64F;

/**
 * @author Marc A. Suchard
 */
public class BranchRateGradient implements GradientWrtParameterProvider, Reportable {

    private final TreeDataLikelihood treeDataLikelihood;
    private final TreeTrait treeTraitProvider;
    private final Tree tree;

    private final int nTaxa;
    private final int nTraits;
    private final int dimTrait;
    private final int dimPartial;

    private final Parameter rateParameter;
    private final MatrixParameterInterface diffusionParameter;
    private final ArbitraryBranchRates branchRateModel;

    public BranchRateGradient(String traitName,
                              TreeDataLikelihood treeDataLikelihood,
                              ContinuousDataLikelihoodDelegate likelihoodDelegate,
                              Parameter rateParameter) {

        assert(treeDataLikelihood != null);

        this.treeDataLikelihood = treeDataLikelihood;
        this.tree = treeDataLikelihood.getTree();
        this.rateParameter = rateParameter;



        this.diffusionParameter = likelihoodDelegate.getDiffusionModel().getPrecisionParameter();
        branchRateModel = (ArbitraryBranchRates) treeDataLikelihood.getBranchRateModel();

        // TODO Move into different constructor / parser
        String fcdName = TipFullConditionalDistributionDelegate.getName(traitName);
        if (treeDataLikelihood.getTreeTrait(fcdName) == null) {
            likelihoodDelegate.addFullConditionalDensityTrait(traitName);
        }

        treeTraitProvider = treeDataLikelihood.getTreeTrait(fcdName);

        assert (treeTraitProvider != null);

        nTaxa = treeDataLikelihood.getTree().getExternalNodeCount();
        nTraits = treeDataLikelihood.getDataLikelihoodDelegate().getTraitCount();
        dimTrait = treeDataLikelihood.getDataLikelihoodDelegate().getTraitDim();
        dimPartial = dimTrait + likelihoodDelegate.getPrecisionType().getMatrixLength(dimTrait);

        if (nTraits != 1) {
            throw new RuntimeException("Not yet implemented for >1 traits");
        }
    }

    @Override
    public Likelihood getLikelihood() {
        return treeDataLikelihood;
    }

    @Override
    public Parameter getParameter() {
        return rateParameter;
    }

    @Override
    public int getDimension() {
        return getParameter().getDimension();
    }

    private class BranchInformation {

        WrappedVector childMean;
        DenseMatrix64F childPrecision;
        DenseMatrix64F childVariance;

        WrappedVector parentMean;
        DenseMatrix64F parentPrecision;
        DenseMatrix64F parentVariance;

        DenseMatrix64F branchPrecision;

        private BranchInformation(final double[] fcd,
                                  final NodeRef node, final NodeRef parent,
                                  final DenseMatrix64F branchPrecision) {

            int nodeOffset = dimPartial * node.getNumber();
            int parentOffset = dimPartial * parent.getNumber();

            this.childMean = new WrappedVector.Raw(fcd, nodeOffset, dimTrait);
            this.childPrecision = MissingOps.wrap(fcd, nodeOffset + dimTrait, dimTrait, dimTrait);
            this.childVariance = MissingOps.wrap(fcd, nodeOffset + dimTrait + dimTrait * dimTrait, dimTrait, dimTrait);

            this.parentMean = new WrappedVector.Raw(fcd, parentOffset, dimTrait);
            this.parentPrecision = MissingOps.wrap(fcd, parentOffset + dimTrait, dimTrait, dimTrait);
            this.parentVariance = MissingOps.wrap(fcd, parentOffset + dimTrait + dimTrait * dimTrait, dimTrait, dimTrait);

            this.branchPrecision = branchPrecision;
        }

        public BranchInformation(double[] fcd, NodeRef node, NodeRef parent) {

        }
    }

    class WrappedMatrix extends DenseMatrix64F {

        public WrappedMatrix() {
            super();
        }
    }

    DenseMatrix64F test;
    
    @Override
    public double[] getGradientLogDensity() {

        double[] result = new double[rateParameter.getDimension()];

        double[] fcd = (double[]) treeTraitProvider.getTrait(tree, null);

        DenseMatrix64F branchPrecision = new DenseMatrix64F(diffusionParameter.getParameterAsMatrix());

        for (int i = 0; i < tree.getNodeCount(); ++i) {
            final NodeRef node = tree.getNode(i);

            if (!tree.isRoot(node)) {

                NodeRef parent = tree.getParent(node);

                BranchInformation branchInformation = new BranchInformation(fcd, node, parent);

                double gradient = 0.0;
                // TODO Compute

                final int destinationIndex = branchRateModel.getParameterIndexFromNode(node);
                assert (destinationIndex != -1);

                result[destinationIndex] = gradient;
            }
        }

        return result;
    }

    private static double getGradientForBranch(BranchInformation branch) {

        double gradient = 0.0;


        if (DEBUG) {
            System.err.println(branch.childMean);
            System.err.println(branch.childPrecision);
            System.err.println(branch.parentMean);
            System.err.println(branch.parentPrecision);
        }

        return gradient;
    }

//    private double[] map(final double[] in) {
//        double[] out = new double[in.length];
//
//        for (int i = 0; i < tree.getNodeCount(); ++i) {
//            NodeRef node = tree.getNode(i);
//            int index = branchRateModel.getParameterIndexFromNode(node);
//            if (index != -1) {
//                out[index] = in[node.getNumber()];
//            }
//        }
//
//        return out;
//    }

    @Override
    public String getReport() {
        return (new dr.math.matrixAlgebra.Vector(getGradientLogDensity())).toString();
    }

    private static final boolean DEBUG = true;
}
