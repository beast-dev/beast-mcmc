/*
 * TreeScaledRepeatedMeasuresTraitDataModel.java
 *
 * Copyright (c) 2002-2019 Alexei Drummond, Andrew Rambaut and Marc Suchard
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
import dr.evomodel.treedatalikelihood.continuous.cdi.PrecisionType;
import dr.inference.model.CompoundParameter;
import dr.inference.model.MatrixParameterInterface;
import org.ejml.data.DenseMatrix64F;
import org.ejml.ops.CommonOps;

/**
 * @author Marc A. Suchard
 * @author Paul Bastide
 */
public class TreeScaledRepeatedMeasuresTraitDataModel extends RepeatedMeasuresTraitDataModel {

    private Tree treeModel;
    private ContinuousRateTransformation rateTransformation;

    public TreeScaledRepeatedMeasuresTraitDataModel(String name,
                                                    ContinuousTraitPartialsProvider childModel,
                                                    CompoundParameter parameter,
                                                    boolean[] missingIndicators,
                                                    boolean useMissingIndices,
                                                    final int dimTrait,
                                                    final int numTraits,
                                                    MatrixParameterInterface samplingPrecision,
                                                    PrecisionType precisionType) {
        super(name, childModel, parameter, missingIndicators, useMissingIndices, dimTrait, numTraits,
                samplingPrecision, precisionType);

        if (!(childModel instanceof ContinuousTraitDataModel)) {
            throw new RuntimeException("not yet implemented for alternative child models. " +
                    "(can't just scale the partial in super.getTipPartial)");
        }
    }

    @Override
    public void addTreeAndRateModel(Tree treeModel, ContinuousRateTransformation rateTransformation) {
        this.treeModel = treeModel;
        this.rateTransformation = rateTransformation;
    }

    @Override
    public double[] getTipPartial(int taxonIndex, boolean fullyObserved) {

        double[] partial = super.getTipPartial(taxonIndex, fullyObserved);
        double tipHeight = getTipHeight(taxonIndex);

        scalePartial(tipHeight, partial);

        if (DEBUG) {
            System.err.println("taxon " + taxonIndex);
            System.err.println("\tscaling: " + tipHeight);
        }

        return partial;
    }

    private double getTipHeight(int taxonIndex) {
        double time = treeModel.getNodeHeight(treeModel.getRoot()) - treeModel.getNodeHeight(treeModel.getExternalNode(taxonIndex));
        return time * rateTransformation.getNormalization();
    }

    private void scalePartial(double tipHeight, double[] partial) {
        scaleArray(1 / tipHeight, partial, dimTrait, dimTrait * dimTrait); // Precision
        scaleArray(tipHeight, partial, dimTrait + dimTrait * dimTrait, dimTrait * dimTrait); // Variance
    }

    private void scaleArray(double t, double[] x, int offset, int length) {
        for (int i = offset; i < offset + length; i++) {
            x[i] *= t;
        }
    }

    @Override
    public void chainRuleWrtVariance(double[] gradient, NodeRef node) {
        double tipHeight = getTipHeight(node.getNumber());
        scaleArray(tipHeight, gradient, 0, gradient.length);
    }

    @Override
    public DenseMatrix64F getExtensionVariance(NodeRef node) {
        DenseMatrix64F var = getExtensionVariance();
        double tipHeight = getTipHeight(node.getNumber());
        CommonOps.scale(tipHeight, var);
        return var;
    }

    @Override
    public void getMeanTipVariances(DenseMatrix64F samplingVariance, DenseMatrix64F samplingComponent) {
        double meanDepth = 0;
        for (int i = 0; i < treeModel.getExternalNodeCount(); i++) {
            meanDepth += getTipHeight(treeModel.getExternalNode(i).getNumber());
        }
        meanDepth /= treeModel.getExternalNodeCount();
        CommonOps.scale(meanDepth, samplingVariance, samplingComponent);
    }

    private static final boolean DEBUG = false;

}
