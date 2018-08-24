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
import dr.evolution.tree.TreeTraitProvider;
import dr.evomodel.treedatalikelihood.*;
import dr.evomodel.treedatalikelihood.preorder.MultivariateConditionalOnTipsRealizedDelegate;
import dr.evomodel.treedatalikelihood.preorder.ProcessSimulationDelegate;
import dr.inference.loggers.LogColumn;
import dr.inference.loggers.Loggable;
import dr.inference.loggers.NumberColumn;
import dr.inference.model.Likelihood;
import dr.inference.model.Parameter;
import dr.math.distributions.GaussianProcessRandomGenerator;

import java.util.List;

import static dr.evomodel.treedatalikelihood.preorder.AbstractRealizedContinuousTraitDelegate.getTipTraitName;
import static dr.evomodel.treedatalikelihood.continuous.ContinuousDataLikelihoodDelegate.createWithMissingData;

/**
 * @author Marc A. Suchard
 */
public class TreeTipGaussianProcess implements GaussianProcessRandomGenerator, Loggable {

    private final String traitName;
    private final ContinuousDataLikelihoodDelegate likelihoodDelegate;
    private final TreeDataLikelihood treeDataLikelihood;
    private final Tree tree;

    private final Parameter maskParameter;

    private final List<Integer> sampleIndices;
    private final int missingLength;
    private final int drawLength;

    private final boolean[] doSample;
    private final boolean[] doNotSample;

    private final TreeTrait tipSampleTrait;

    private final boolean truncateToMissingOnly;

    private final LikelihoodTreeTraversal treeTraversalDelegate;

    public TreeTipGaussianProcess(String traitName,
                                  TreeDataLikelihood treeDataLikelihood,
                                  ContinuousDataLikelihoodDelegate likelihoodDelegate,
                                  Parameter maskParameter,
                                  boolean truncateToMissingOnly) {

        this.treeDataLikelihood = treeDataLikelihood;
        this.tree = treeDataLikelihood.getTree();
        this.maskParameter = maskParameter;

        List<Integer> indices = likelihoodDelegate.getDataModel().getMissingIndices();
        if (indices == null || indices.size() == 0) {
            traitName = traitName + ".missing";
            likelihoodDelegate = createWithMissingData(likelihoodDelegate);

            ProcessSimulationDelegate simulationDelegate =
                    new MultivariateConditionalOnTipsRealizedDelegate(traitName, treeDataLikelihood.getTree(),
                            likelihoodDelegate.getDiffusionModel(), likelihoodDelegate.getDataModel(), likelihoodDelegate.getRootPrior(),
                            likelihoodDelegate.getRateTransformation(), likelihoodDelegate);

            TreeTraitProvider traitProvider = new ProcessSimulation(treeDataLikelihood, simulationDelegate);

            treeDataLikelihood.addTraits(traitProvider.getTreeTraits());

            this.treeTraversalDelegate = new LikelihoodTreeTraversal(
                    treeDataLikelihood.getTree(),
                    treeDataLikelihood.getBranchRateModel(),
                    TreeTraversal.TraversalType.POST_ORDER);

        } else {
            this.treeTraversalDelegate = null; // No extra post-order traversal necessary
        }

        this.likelihoodDelegate = likelihoodDelegate;
        this.sampleIndices = likelihoodDelegate.getDataModel().getMissingIndices();
        this.missingLength = sampleIndices.size();
        this.traitName = traitName;

        String tipTraitName = getTipTraitName(traitName);
        tipSampleTrait = treeDataLikelihood.getTreeTrait(tipTraitName);

        assert (tipSampleTrait != null);

        double[] draw = drawAllTraits();
        this.drawLength = draw.length;

        doSample = new boolean[drawLength];
        for (int i : sampleIndices) {
            doSample[i] = true;
        }

        doNotSample = new boolean[drawLength];
        for (int i = 0; i < drawLength; ++i) {
            doNotSample[i] = !doSample[i];
        }

        this.truncateToMissingOnly = truncateToMissingOnly;
    }

    private double[] drawAllTraits() {
        //treeDataLikelihood.fireModelChanged();
        doPostOrderTraversal();
        return (double []) tipSampleTrait.getTrait(treeDataLikelihood.getTree(), null);
    }

    private double[] drawSelectedTraits() {
        double[] draw = drawAllTraits();
        if (truncateToMissingOnly) {
             draw = maskDraw(draw, doSample);
         }
         return draw;
    }

    private void doPostOrderTraversal() {

        if (treeTraversalDelegate != null) {
            
            treeTraversalDelegate.updateAllNodes();
            treeTraversalDelegate.dispatchTreeTraversalCollectBranchAndNodeOperations();

            final List<DataLikelihoodDelegate.BranchOperation> branchOperations = treeTraversalDelegate.getBranchOperations();
            final List<DataLikelihoodDelegate.NodeOperation> nodeOperations = treeTraversalDelegate.getNodeOperations();

            final NodeRef root = treeDataLikelihood.getTree().getRoot();

            try {
                likelihoodDelegate.calculateLikelihood(branchOperations, nodeOperations, root.getNumber());

            } catch (DataLikelihoodDelegate.LikelihoodException e) {
                throw new RuntimeException("Unhandled exception");
            }
        } else {
            treeDataLikelihood.fireModelChanged();
        }
    }

//    private double[] crop(double[] in, int length) {
//        double[] out = new double[length];
//        System.arraycopy(in, 0, out, 0, length);
//        return out;
//    }

    private double[] maskDraw(double[] in, boolean[] mask) {
        double[] out = new double[missingLength];

        int index = 0;
        for (int i = 0; i < in.length; ++i) {
            if (mask[i]) {
                out[index] = in[i];
                ++index;
            }
        }

        return out;
    }

    @Override
    public Likelihood getLikelihood() {
        return treeDataLikelihood;
    }

    @Override
    public int getDimension() {
        return truncateToMissingOnly ? missingLength : drawLength;
    }

    @Override
    public double[][] getPrecisionMatrix() {
        throw new RuntimeException("Precision matrix is never formed");
    }

    @Override
    public Object nextRandom() {
        return drawSelectedTraits();
    }

    @Override
    public double logPdf(Object x) {
        throw new RuntimeException("Density is never evaluated");
    }

    @Override
    public LogColumn[] getColumns() {

        double[] sample = drawSelectedTraits();

        LogColumn[] columns = new LogColumn[sample.length];

        for (int i = 0; i < sample.length; ++i) {
            columns[i] = new TipSampleColumn("tip", i);
        }

        return columns;
    }

    private double[] currentSample = null;

    private class TipSampleColumn extends NumberColumn {

        private int index;

        TipSampleColumn(String label, int index) {
            super(label);
            this.index = index;
        }

        @Override
        public double getDoubleValue() {

            if (index == 0) {
                currentSample = drawSelectedTraits();
            }

            return currentSample[index];
        }
    }
}
