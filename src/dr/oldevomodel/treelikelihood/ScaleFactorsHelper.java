/*
 * ScaleFactorsHelper.java
 *
 * Copyright (c) 2002-2015 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.oldevomodel.treelikelihood;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evomodel.treelikelihood.LikelihoodScalingProvider;
import dr.inference.model.AbstractModel;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;

import java.util.Arrays;

/**
 * @author Marc A. Suchard
 */
@Deprecated // Switching to BEAGLE
public class ScaleFactorsHelper extends AbstractModel {

    public ScaleFactorsHelper(LikelihoodScalingProvider scalingProvider,
                              Model treeLikelihood,
                              Tree treeModel,
                              int stateCount,
                              int patternCount,
                              int categoryCount) {
        super("ScaleFactorHelper");
        this.scalingProvider = scalingProvider;
        this.treeLikelihood = treeLikelihood;
        this.treeModel = treeModel;
        this.stateCount = stateCount;
        this.patternCount = patternCount;
        this.categoryCount = categoryCount;

        addModel(treeLikelihood);
    }

    private double[][] scaleFactors = null;
    private double[][] storedScaleFactors = null;
    private boolean scaleFactorsKnown = false;
    private boolean storedScaleFactorsKnown = false;
    private double[] buffer = null;

    final private LikelihoodScalingProvider scalingProvider;
    final private Model treeLikelihood;
    final private Tree treeModel;
    final private int stateCount;
    final private int patternCount;
    final private int categoryCount;

    private void addScaleFactors(final double[] in0, final double[] in1, final double[] in2,
                                 double[] out, final int length) {
        for (int i = 0; i < length; ++i) {
            out[i] = in0[i] + in1[i] + in2[i];
        }
    }

    private void traverseComputeScaleFactors(final Tree tree, final NodeRef node) {

        final int nodeNumber = node.getNumber();

        if (tree.isExternal(node)) {
            Arrays.fill(scaleFactors[nodeNumber], 0.0);
        } else {
            final NodeRef child0 = tree.getChild(node, 0);
            final NodeRef child1 = tree.getChild(node, 1);

            traverseComputeScaleFactors(tree, child0);
            traverseComputeScaleFactors(tree, child1);

            if (!tree.isExternal(child0) || !tree.isExternal(child1)) {
                scalingProvider.getLogScalingFactors(nodeNumber, buffer);
                addScaleFactors(scaleFactors[child0.getNumber()], scaleFactors[child1.getNumber()],
                        buffer, scaleFactors[nodeNumber], patternCount);
            } else {
                scalingProvider.getLogScalingFactors(nodeNumber, scaleFactors[nodeNumber]);
            }
        }
    }

    private void computeScaleFactors() {
        if (scaleFactors == null) {
            scaleFactors = new double[treeModel.getNodeCount()][patternCount];
            storedScaleFactors = new double[treeModel.getNodeCount()][patternCount];
        }
        if (buffer == null) {
            buffer = new double[patternCount];
        }
        traverseComputeScaleFactors(treeModel, treeModel.getRoot());
    }

    public void resetScaleFactors() {
        scaleFactorsKnown = false;
    }

    public void rescalePartials(final int nodeNumber, double[] partials) {
        if (scalingProvider.arePartialsRescaled() && nodeNumber >= treeModel.getExternalNodeCount()) {

            if (!scaleFactorsKnown) {
                computeScaleFactors();
                scaleFactorsKnown = true;
            }

            int index = 0;
            for (int category = 0; category < categoryCount; ++category) {
                for (int pattern = 0; pattern < patternCount; ++pattern) {
                    final double scale = Math.exp(scaleFactors[nodeNumber][pattern]);
                    for (int state = 0; state < stateCount; ++state) {
                        partials[index] *= scale;
                        index++;
                    }
                }
            }
        }
    }

    protected void handleModelChangedEvent(Model model, Object object, int index) {
//        if (model == treeLikelihood) {
//            scaleFactorsKnown = false; // TODO Not all tree likelihood changes affect scale factors
//        } else {
//            throw new IllegalArgumentException("Illegal model argument");
//        }
    }

    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        // Do nothing
    }

    protected void storeState() {
//        for (int i = 0; i < scaleFactors.length; ++i) {
//            System.arraycopy(scaleFactors[i], 0, storedScaleFactors[i], 0, scaleFactors[i].length);
//        }
//        storedScaleFactorsKnown = scaleFactorsKnown;
    }

    protected void restoreState() {
//        double[][] tmp = storedScaleFactors;
//        storedScaleFactors = scaleFactors;
//        scaleFactors = tmp;
//
//        scaleFactorsKnown = storedScaleFactorsKnown;
    }

    protected void acceptState() {
        // Do nothing
    }
}

