/*
 * BastaLikelihoodDelegate.java
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

import dr.evolution.coalescent.IntervalType;
import dr.evolution.tree.Tree;
import dr.evomodel.bigfasttree.BigFastTreeIntervals;
import dr.evomodel.substmodel.EigenDecomposition;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.*;
import dr.math.matrixAlgebra.WrappedVector;
import dr.util.Citable;
import dr.util.Citation;
import dr.xml.Reportable;

import java.util.ArrayList;
import java.util.List;

/**
 * BastaLikelihoodDelegate - interface for a plugin delegate for the BASTA model likelihood.
 *
 * @author Marc A. Suchard
 * @author Yucai Shao
 * @author Guy Baele
 */
public interface BastaLikelihoodDelegate extends ProcessOnCoalescentIntervalDelegate, Model, Profileable, Reportable {

    void makeDirty();

    void storeState();

    void restoreState();

    double calculateLikelihood(List<BranchIntervalOperation> branchOperations,
                               List<TransitionMatrixOperation> matrixOperations,
                               List<Integer> intervalStarts,
                               int rootNodeNumber);

    default void setPartials(int index, double[] partials) {
        throw new RuntimeException("Not yet implemented");
    }

    default void getPartials(int index, double[] partials) {
        assert index >= 0;
        assert partials != null;

        throw new RuntimeException("Not yet implemented");
    }

    default void updateEigenDecomposition(int index, EigenDecomposition decomposition, boolean flip) {
        throw new RuntimeException("Not yet implemented");
    }

    default void updatePopulationSizes(int index, double[] sizes, boolean flip) {
        throw new RuntimeException("Not yet implemented");
    }

    double[][] calculateGradient(List<BranchIntervalOperation> branchOperations,
                                 List<TransitionMatrixOperation> matrixOperation,
                                 List<Integer> intervalStarts,
                                 int rootNodeNumber);

    double[] calculateGradientPopSize(List<BranchIntervalOperation> branchOperations, List<TransitionMatrixOperation> matrixOperations, List<Integer> intervalStarts, int number);

    abstract class AbstractBastaLikelihoodDelegate extends AbstractModel implements BastaLikelihoodDelegate, Citable {

        protected static final boolean PRINT_COMMANDS = false;

        protected final int maxNumCoalescentIntervals;

        protected final ParallelizationScheme parallelizationScheme;

        protected final int stateCount;

        protected final Tree tree;

        protected final boolean transpose;

        public AbstractBastaLikelihoodDelegate(String name,
                                               Tree tree,
                                               int stateCount, boolean transpose) {
            super(name);

            this.tree = tree;
            this.stateCount = stateCount;
            this.maxNumCoalescentIntervals = getMaxNumberOfCoalescentIntervals(tree);
            this.transpose = transpose;
            this.parallelizationScheme = ParallelizationScheme.NONE;
        }

        private int getMaxNumberOfCoalescentIntervals(Tree tree) {
            BigFastTreeIntervals intervals = new BigFastTreeIntervals((TreeModel) tree); // TODO fix BFTI to take a Tree
            int zeroLengthSampling = 0;
            for (int i = 0; i < intervals.getIntervalCount(); ++i) {
                if (intervals.getIntervalType(i) == IntervalType.SAMPLE && intervals.getIntervalTime(i) == 0.0) {
                    ++zeroLengthSampling;
                }
            }
            return tree.getNodeCount() - zeroLengthSampling;
        }

        @Override
        public void makeDirty() {

        }

        @Override
        protected void handleModelChangedEvent(Model model, Object object, int index) {
            throw new RuntimeException("Should not be called");
        }

        @Override
        protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
            throw new RuntimeException("Should not be called");
        }

        @Override
        public void storeState() {

        }

        @Override
        public void restoreState() {

        }

        @Override
        protected void acceptState() {

        }

        enum ParallelizationScheme {
            NONE,
            FULL
        }

        abstract protected void computeBranchIntervalOperations(List<Integer> intervalStarts,
                                                                List<BranchIntervalOperation> branchIntervalOperations);

        abstract protected void computeTransitionProbabilityOperations(List<TransitionMatrixOperation> matrixOperations);

        abstract protected double computeCoalescentIntervalReduction(List<Integer> intervalStarts,
                                                                     List<BranchIntervalOperation> branchIntervalOperations);

        @Override
        public double calculateLikelihood(List<BranchIntervalOperation> branchOperations,
                                          List<TransitionMatrixOperation> matrixOperation,
                                          List<Integer> intervalStarts,
                                          int rootNodeNumber) {

            if (PRINT_COMMANDS) {
                System.err.println("Tree = " + tree);
            }

            computeTransitionProbabilityOperations(matrixOperation);
            computeBranchIntervalOperations(intervalStarts, branchOperations);

            double logL = computeCoalescentIntervalReduction(intervalStarts, branchOperations);

            if (PRINT_COMMANDS) {
                System.err.println("logL = " + logL + " " + getStamp() + "\n");
                if (printCount > 1000) {
                    System.exit(-1);
                }
                ++printCount;
            }

            return logL;
        }

        abstract protected void computeBranchIntervalOperationsGrad(List<Integer> intervalStarts, List<TransitionMatrixOperation> matrixOperations,
                                                                List<BranchIntervalOperation> branchIntervalOperations);

        abstract protected void computeTransitionProbabilityOperationsGrad(List<TransitionMatrixOperation> matrixOperations);

        abstract protected double[][] computeCoalescentIntervalReductionGrad(List<Integer> intervalStarts,
                                                                     List<BranchIntervalOperation> branchIntervalOperations);

        abstract protected double[] computeCoalescentIntervalReductionGradPopSize(List<Integer> intervalStarts,
                                                                             List<BranchIntervalOperation> branchIntervalOperations);

        @Override
        public double[][] calculateGradient(List<BranchIntervalOperation> branchOperations,
                                            List<TransitionMatrixOperation> matrixOperation,
                                            List<Integer> intervalStarts,
                                            int rootNodeNumber) {

            if (PRINT_COMMANDS) {
                System.err.println("Tree = " + tree);
            }

            computeTransitionProbabilityOperationsGrad(matrixOperation);

            computeBranchIntervalOperationsGrad(intervalStarts, matrixOperation, branchOperations);

            double[][] grad = computeCoalescentIntervalReductionGrad(intervalStarts, branchOperations);

//            if (PRINT_COMMANDS) {
//                System.err.println("logL = " + logL + " " + getStamp() + "\n");
//                if (printCount > 1000) {
//                    System.exit(-1);
//                }
//                ++printCount;
//            }

            return grad;
        }

        @Override
        public double[]calculateGradientPopSize(List<BranchIntervalOperation> branchOperations,
                                            List<TransitionMatrixOperation> matrixOperation,
                                            List<Integer> intervalStarts,
                                            int rootNodeNumber) {

            if (PRINT_COMMANDS) {
                System.err.println("Tree = " + tree);
            }

            computeBranchIntervalOperationsGrad(intervalStarts, matrixOperation, branchOperations);

            double[] grad = computeCoalescentIntervalReductionGradPopSize(intervalStarts, branchOperations);

//            if (PRINT_COMMANDS) {
//                System.err.println("logL = " + logL + " " + getStamp() + "\n");
//                if (printCount > 1000) {
//                    System.exit(-1);
//                }
//                ++printCount;
//            }

            return grad;
        }


        abstract String getStamp();

        int printCount = 0;

        @Override
        public long getTotalCalculationCount() {
            return 0;
        }

        @Override
        public Citation.Category getCategory() {
            return null;
        }

        @Override
        public String getDescription() {
            return null;
        }

        @Override
        public List<Citation> getCitations() {
            return new ArrayList<>();
        }

        @Override
        public String getReport() {
            return null;
        }
    }
}
