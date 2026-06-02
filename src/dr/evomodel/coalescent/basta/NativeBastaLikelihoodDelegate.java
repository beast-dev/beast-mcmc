/*
 * NativeBastaLikelihoodDelegate.java
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

package dr.evomodel.coalescent.basta;

import dr.evolution.alignment.PatternList;
import dr.evolution.tree.Tree;
import dr.evomodel.substmodel.SubstitutionModel;
import dr.inference.model.Parameter;

import java.util.List;

/**
 * @author Marc A. Suchard
 */
public class NativeBastaLikelihoodDelegate extends BastaLikelihoodDelegate.AbstractBastaLikelihoodDelegate {

    private final NativeBastaJniWrapper jni;

    public NativeBastaLikelihoodDelegate(String name,
                                         Tree tree,
                                         int stateCount,
                                         boolean transpose) {
        super(name, tree, stateCount, transpose);
        jni = NativeBastaJniWrapper.getBastaJniWrapper();
    }

    @Override
    String getStamp() { return "native"; }


    @Override
    protected void computeBranchIntervalOperations(List<Integer> intervalStarts,
                                                   List<BranchIntervalOperation> branchIntervalOperations) {
        if (branchIntervalOperations != null) {
            for (BranchIntervalOperation operation : branchIntervalOperations) {
                System.err.println(operation.toString());
            }
        }
    }

    @Override
    protected void computeTransitionProbabilityOperations(List<TransitionMatrixOperation> matrixOperations) {
        if (matrixOperations != null) {
            for (TransitionMatrixOperation operation : matrixOperations) {
                System.err.println(operation.toString());
            }
        }
    }

    @Override
    protected double computeCoalescentIntervalReduction(List<Integer> intervalStarts,
                                                        List<BranchIntervalOperation> branchIntervalOperations) {
        if (intervalStarts != null) {
            for (int start : intervalStarts) {
                System.err.println(start);
            }
        }

        return 0.0;
    }

    @Override
    protected void computeBranchIntervalOperationsGrad(List<Integer> intervalStarts, List<TransitionMatrixOperation> matrixOperations, List<BranchIntervalOperation> branchIntervalOperations) {

    }

    @Override
    protected void computeTransitionProbabilityOperationsGrad(List<TransitionMatrixOperation> matrixOperations) {

    }

    @Override
    protected double[][] computeCoalescentIntervalReductionGrad(List<Integer> intervalStarts, List<BranchIntervalOperation> branchIntervalOperations) {
        return new double[0][];
    }

    @Override
    protected double[] computeCoalescentIntervalReductionGradPopSize(List<Integer> intervalStarts, List<BranchIntervalOperation> branchIntervalOperations) {
        return new double[0];
    }
}
