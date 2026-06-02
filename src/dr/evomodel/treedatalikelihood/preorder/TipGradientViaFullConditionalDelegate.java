/*
 * TipGradientViaFullConditionalDelegate.java
 *
 * Copyright © 2002-2024 the BEAST Development Team
 * http://beast.community/about
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
 *
 */

package dr.evomodel.treedatalikelihood.preorder;

import dr.evolution.tree.MutableTreeModel;
import dr.evolution.tree.NodeRef;
import dr.evomodel.continuous.MultivariateDiffusionModel;
import dr.evomodel.treedatalikelihood.continuous.*;
import dr.evomodel.treedatalikelihood.continuous.cdi.PrecisionType;
import dr.inference.model.MatrixParameterInterface;
import dr.math.matrixAlgebra.missingData.MissingOps;
import org.ejml.data.DenseMatrix64F;
import org.ejml.ops.CommonOps;

import static dr.math.matrixAlgebra.missingData.MissingOps.safeInvert2;

/**
 * @author Marc A. Suchard
 */
public class TipGradientViaFullConditionalDelegate extends TipFullConditionalDistributionDelegate {

    private final int offset;
    private final int dimGradient;

    private final int[] subInds;

    private final boolean doSubset;

    public TipGradientViaFullConditionalDelegate(String name, MutableTreeModel tree,
                                                 MultivariateDiffusionModel diffusionModel,
                                                 ContinuousTraitPartialsProvider dataModel,
                                                 ConjugateRootTraitPrior rootPrior,
                                                 ContinuousRateTransformation rateTransformation,
                                                 ContinuousDataLikelihoodDelegate likelihoodDelegate,
                                                 int offset,
                                                 int dimGradient) {

        super(name, tree, diffusionModel, dataModel, rootPrior, rateTransformation, likelihoodDelegate);
        this.offset = offset;
        this.dimGradient = dimGradient;

        subInds = new int[dimGradient];
        for (int i = 0; i < dimGradient; i++) {
            subInds[i] = i + offset;
        }

        if (offset != 0 || dimGradient != dimTrait) {
            doSubset = true;
        } else {
            doSubset = false;
        }
    }

    public static String getName(String name) {
        return "grad." + name;
    }

    public String getTraitName(String name) {
        return getName(name);
    }

    @Override
    protected double[] getTraitForNode(NodeRef node) {

        if (likelihoodDelegate.getPrecisionType() == PrecisionType.SCALAR) {
            return getTraitForNodeScalar(node);
        } else if (likelihoodDelegate.getPrecisionType() == PrecisionType.FULL) {
            return getTraitForNodeFull(node);
        } else {
            throw new RuntimeException("Tip gradients are not implemented for '" +
                    likelihoodDelegate.getPrecisionType().toString() + "' likelihoods");
        }
    }

    private double[] getTraitForNodeScalar(NodeRef node) {

        final double[] fullConditionalPartial = super.getTraitForNode(node);

        final double[] postOrderPartial = new double[dimPartial * numTraits];
        cdi.getPostOrderPartial(likelihoodDelegate.getActiveNodeIndex(node.getNumber()), postOrderPartial);

        final MatrixParameterInterface precision = diffusionModel.getPrecisionParameter();

        final double[] gradient = new double[dimTrait * numTraits];

        if (numTraits > 1) {
            throw new RuntimeException("Not yet implemented");
        }

        final double scale = fullConditionalPartial[dimTrait];

        for (int i = 0; i < dimTrait; ++i) {

            double sum = 0.0;
            for (int j = 0; j < dimTrait; ++j) {
                sum += (fullConditionalPartial[j] - postOrderPartial[j]) * scale *
                        precision.getParameterValue(i, j);
//                        precision.getParameterValue(i * dimTrait + j);
            }

            gradient[i] = sum;
        }

        return gradient;
    }

    protected double[] getTraitForNodeFull(NodeRef node) {

        if (numTraits > 1) {
            throw new RuntimeException("Not yet implemented");
        }

        // Pre stats
        final double[] fullConditionalPartial = super.getTraitForNode(node);
        NormalSufficientStatistics statPre = new NormalSufficientStatistics(fullConditionalPartial, 0, dimTrait, Pd, likelihoodDelegate.getPrecisionType());
        DenseMatrix64F precisionPre = statPre.getRawPrecisionCopy();
        DenseMatrix64F meanPre = statPre.getRawMeanCopy();

        // Post mean
        final double[] postOrderPartial = new double[dimPartial * numTraits];
        int nodeIndex = likelihoodDelegate.getActiveNodeIndex(node.getNumber());
        cdi.getPostOrderPartial(nodeIndex, postOrderPartial);
        DenseMatrix64F meanPost = MissingOps.wrap(postOrderPartial, 0, dimTrait, 1);

        if (doSubset) {
            precisionPre = MissingOps.gatherRowsAndColumns(precisionPre, subInds, subInds);
            meanPre = MissingOps.gatherRowsAndColumns(meanPre, subInds, new int[]{0});
            meanPost = MissingOps.gatherRowsAndColumns(meanPost, subInds, new int[]{0});
        }

        if (dataModel.needToUpdateTipDataGradient(offset, dimGradient)) {
            DenseMatrix64F variancePre = statPre.getRawVarianceCopy();
            if (doSubset) {
                variancePre = MissingOps.gatherRowsAndColumns(precisionPre, subInds, subInds);
            }

            dataModel.updateTipDataGradient(precisionPre, variancePre, node, offset, dimGradient);
        }

        // - Q_i * (X_i - m_i)
        DenseMatrix64F gradient = new DenseMatrix64F(dimGradient, numTraits);
        CommonOps.addEquals(meanPost, -1.0, meanPre);
        CommonOps.changeSign(meanPost);
        CommonOps.mult(precisionPre, meanPost, gradient);

        return gradient.getData();
    }
}
