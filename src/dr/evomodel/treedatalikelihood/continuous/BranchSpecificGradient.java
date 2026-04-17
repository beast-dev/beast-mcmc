/*
 * BranchSpecificGradient.java
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

package dr.evomodel.treedatalikelihood.continuous;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeTrait;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.evomodel.treedatalikelihood.preorder.BranchConditionalDistributionDelegate;
import dr.evomodel.treedatalikelihood.preorder.BranchSufficientStatistics;
import dr.evomodel.treedatalikelihood.preorder.NormalSufficientStatistics;
import dr.inference.hmc.GradientWrtParameterProvider;
import dr.inference.loggers.LogColumn;
import dr.inference.loggers.Loggable;
import dr.inference.model.Likelihood;
import dr.inference.model.Parameter;
import dr.math.MultivariateFunction;
import dr.math.NumericalDerivative;
import dr.xml.Reportable;
import org.ejml.data.DenseMatrix64F;
import org.ejml.ops.CommonOps;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * @author Marc A. Suchard
 * @author Paul Bastide
 */
public class BranchSpecificGradient implements GradientWrtParameterProvider, Reportable, Loggable {

    private final TreeDataLikelihood treeDataLikelihood;
    private final TreeTrait<List<BranchSufficientStatistics>> treeTraitProvider;
    private final Tree tree;
    private final int nTraits;
    //    private final int dim;
    private final Parameter parameter;
    //    private final ArbitraryBranchRates branchRateModel;
    private final ContinuousTraitGradientForBranch branchProvider;
    private static final String BRANCH_DEBUG_PROPERTY = "beast.debug.branchGradient";
    private static final String BRANCH_DEBUG_ONLY_PARAMETER_PROPERTY = "beast.debug.branchGradient.onlyParameterContains";
    private static final String BRANCH_DEBUG_NODE_FILTER_PROPERTY = "beast.debug.branchGradient.nodeFilter";
    private static final String ENABLE_GLOBAL_REMAINDER_PASS_PROPERTY =
            "beast.experimental.enableGlobalRemainderPass";

    public BranchSpecificGradient(String traitName,
                                  TreeDataLikelihood treeDataLikelihood,
                                  ContinuousDataLikelihoodDelegate likelihoodDelegate,
                                  ContinuousTraitGradientForBranch branchProvider,
                                  Parameter parameter) {

        assert (treeDataLikelihood != null);

        this.treeDataLikelihood = treeDataLikelihood;
        this.tree = treeDataLikelihood.getTree();
        this.parameter = parameter;
        this.branchProvider = branchProvider;

        if (likelihoodDelegate.usesCanonicalOULikelihood()) {
            throw new IllegalArgumentException(
                    "BranchSpecificGradient is not available when traitDataLikelihood uses the canonical OU backend. "
                            + "Use a canonical XML gradient implementation instead.");
        }

        // TODO Move into different constructor / parser
        String bcdName = BranchConditionalDistributionDelegate.getName(traitName);
        if (treeDataLikelihood.getTreeTrait(bcdName) == null) {
            likelihoodDelegate.addBranchConditionalDensityTrait(traitName);
        }

        @SuppressWarnings("unchecked")
        TreeTrait<List<BranchSufficientStatistics>> unchecked = treeDataLikelihood.getTreeTrait(bcdName);
        treeTraitProvider = unchecked;

        assert (treeTraitProvider != null);

        nTraits = treeDataLikelihood.getDataLikelihoodDelegate().getTraitCount();
        if (nTraits != 1) {
            throw new RuntimeException("Not yet implemented for >1 traits");
        }

//        dim = treeDataLikelihood.getDataLikelihoodDelegate().getTraitDim();

    }

    @Override
    public Likelihood getLikelihood() {
        return treeDataLikelihood;
    }

    @Override
    public Parameter getParameter() {
        return parameter;
    }

    @Override
    public int getDimension() {
        return getParameter().getDimension();
    }

    @Override
    public double[] getGradientLogDensity() {

        int dimGradient = branchProvider.getDimension();
        // Ensure branch sufficient statistics are evaluated at the current model state.
        // This avoids one-call stale gradients after parameter updates during HMC checks.
        treeDataLikelihood.getLogLikelihood();

        double[] result = new double[parameter.getDimension()];
        final boolean debugPerBranchComparison =
                shouldEmitBranchDebugForCurrentParameter(Boolean.getBoolean(BRANCH_DEBUG_PROPERTY));
        final double[][] perNodeGradient = debugPerBranchComparison ? new double[tree.getNodeCount()][] : null;
        final double[][] perNodeNumericalGradient = debugPerBranchComparison ? new double[tree.getNodeCount()][] : null;
        final int[] destinationIndexByNode = debugPerBranchComparison ? new int[tree.getNodeCount()] : null;

        // TODO Do single call to traitProvider with node == null (get full tree)
//        List<BranchSufficientStatistics> statisticsForTree = (List<BranchSufficientStatistics>)
//                treeTraitProvider.getTrait(tree, null);

        for (int i = 0; i < tree.getNodeCount(); ++i) {
            final NodeRef node = tree.getNode(i);
            if (!shouldIncludeNodeForAccumulation(node)) {
                continue;
            }

//            if (!tree.isRoot(node)) {

            List<BranchSufficientStatistics> statisticsForNode = treeTraitProvider.getTrait(tree, node);

            assert (statisticsForNode.size() == nTraits);

            double[] gradient;
//                for (int trait = 0; trait < nTraits; ++trait) { // TODO deal with several traits
            try {
                gradient = branchProvider.getGradientForBranch(statisticsForNode.get(0), node);
            } catch (RuntimeException ex) {
                final BranchSufficientStatistics failingStats = statisticsForNode.get(0);
                throw new IllegalStateException(
                        "Branch gradient failure"
                                + " node=" + node.getNumber()
                                + " isRoot=" + tree.isRoot(node)
                                + " isExternal=" + tree.isExternal(node)
                                + " parameterName=" + parameter.getParameterName()
                                + " parameterValues=" + Arrays.toString(parameter.getParameterValues())
                                + " statsSummary={" + summarizeBranchStatistics(failingStats) + "}",
                        ex);
            }
//                }

            final int destinationIndex = getParameterIndexFromNode(node);
            assert (destinationIndex != -1);
            accumulateBranchGradient(result, gradient, destinationIndex, dimGradient);
            if (debugPerBranchComparison) {
                perNodeGradient[node.getNumber()] = Arrays.copyOf(gradient, gradient.length);
                final double[] numericalGradientForNode = branchProvider.getNumericalGradientForBranch(statisticsForNode.get(0), node);
                if (numericalGradientForNode != null) {
                    perNodeNumericalGradient[node.getNumber()] =
                            Arrays.copyOf(numericalGradientForNode, numericalGradientForNode.length);
                }
                destinationIndexByNode[node.getNumber()] = destinationIndex;
            }
//            }
        }

        if (Boolean.getBoolean(ENABLE_GLOBAL_REMAINDER_PASS_PROPERTY)) {
            // Second pass: add global pruning-remainder contributions.
            // For non-selection providers this is a no-op.
            for (int i = 0; i < tree.getNodeCount(); ++i) {
                final NodeRef node = tree.getNode(i);
                if (tree.isRoot(node)) {
                    continue;
                }
                final NodeRef parent = tree.getParent(node);
                final List<BranchSufficientStatistics> statisticsForNode = treeTraitProvider.getTrait(tree, node);
                final BranchSufficientStatistics branchStatistics = statisticsForNode.get(0);
                final double[] globalGradient = new double[dimGradient];
                branchProvider.accumulateGlobalRemainderGradientForBranch(
                        branchStatistics, node, parent, globalGradient);
                final int destinationIndex = getParameterIndexFromNode(node);
                assert (destinationIndex != -1);
                accumulateBranchGradient(result, globalGradient, destinationIndex, dimGradient);
            }
        }

        if (debugPerBranchComparison) {
            emitPerBranchGradientDebug(
                    result,
                    perNodeGradient,
                    perNodeNumericalGradient,
                    destinationIndexByNode,
                    dimGradient);
        }

        return result;
    }

    private static String summarizeBranchStatistics(final BranchSufficientStatistics statistics) {
        final StringBuilder sb = new StringBuilder(512);
        sb.append("branchActualization=").append(finitenessSummary(statistics.getBranch().getRawActualization()));
        sb.append(", branchDisplacement=").append(finitenessSummary(statistics.getBranch().getRawDisplacement()));
        sb.append(", branchPrecision=").append(finitenessSummary(statistics.getBranch().getRawPrecision()));
        sb.append(", aboveMean=").append(finitenessSummary(statistics.getAbove().getRawMean()));
        sb.append(", abovePrecision=").append(finitenessSummary(statistics.getAbove().getRawPrecision()));
        sb.append(", belowMean=").append(finitenessSummary(statistics.getBelow().getRawMean()));
        sb.append(", belowPrecision=").append(finitenessSummary(statistics.getBelow().getRawPrecision()));
        sb.append(", missing=").append(Arrays.toString(statistics.getMissing()));
        return sb.toString();
    }

    private static String finitenessSummary(final DenseMatrix64F matrix) {
        if (matrix == null) {
            return "null";
        }
        int nan = 0;
        int inf = 0;
        double min = Double.POSITIVE_INFINITY;
        double max = Double.NEGATIVE_INFINITY;
        final double[] data = matrix.getData();
        for (double value : data) {
            if (Double.isNaN(value)) {
                nan++;
            } else if (Double.isInfinite(value)) {
                inf++;
            } else {
                if (value < min) {
                    min = value;
                }
                if (value > max) {
                    max = value;
                }
            }
        }
        if (nan + inf == data.length) {
            min = Double.NaN;
            max = Double.NaN;
        }
        return "rows=" + matrix.numRows
                + ",cols=" + matrix.numCols
                + ",nan=" + nan
                + ",inf=" + inf
                + ",min=" + min
                + ",max=" + max;
    }

    private boolean shouldIncludeNodeForAccumulation(final NodeRef node) {
        final String filter = System.getProperty(BRANCH_DEBUG_NODE_FILTER_PROPERTY, "all");
        if ("all".equalsIgnoreCase(filter)) {
            return true;
        }
        if ("rootOnly".equalsIgnoreCase(filter)) {
            return tree.isRoot(node);
        }
        if ("nonRootOnly".equalsIgnoreCase(filter)) {
            return !tree.isRoot(node);
        }
        throw new IllegalArgumentException(
                "Unknown value for " + BRANCH_DEBUG_NODE_FILTER_PROPERTY + ": " + filter +
                        " (supported: all, rootOnly, nonRootOnly)");
    }

    /**
     * Computes the branch-specific gradient after explicitly synchronizing the
     * likelihood state at the current parameter values.
     *
     * This is an opt-in hook for debugging/experimental gradient providers that
     * need a fresh likelihood snapshot before branch-statistics assembly.
     */
    public double[] getGradientLogDensityFromFreshLikelihoodSnapshot() {
        treeDataLikelihood.makeDirty();
        treeDataLikelihood.getLogLikelihood();
        return getGradientLogDensity();
    }

    private void accumulateBranchGradient(double[] destination,
                                          double[] localGradient,
                                          int destinationIndex,
                                          int blockDimension) {
        if (localGradient.length != blockDimension) {
            throw new IllegalStateException(
                    "Unexpected branch local gradient length: expected " + blockDimension +
                            " but found " + localGradient.length + "."
            );
        }
        for (int j = 0; j < blockDimension; j++) {
            final int resultIndex = getResultIndex(destinationIndex, j, blockDimension, destination.length);
            destination[resultIndex] += localGradient[j];
        }
    }

    private int getResultIndex(int destinationIndex, int localIndex, int blockDimension, int resultLength) {
        final int resultIndex = destinationIndex * blockDimension + localIndex;
        if (resultIndex < 0 || resultIndex >= resultLength) {
            throw new IllegalStateException(
                    "Branch gradient index out of bounds: destinationIndex=" + destinationIndex +
                            ", localIndex=" + localIndex +
                            ", blockDimension=" + blockDimension +
                            ", resultLength=" + resultLength + "."
            );
        }
        return resultIndex;
    }

    private boolean shouldEmitBranchDebugForCurrentParameter(final boolean debugEnabled) {
        if (!debugEnabled) {
            return false;
        }
        final String filter = System.getProperty(BRANCH_DEBUG_ONLY_PARAMETER_PROPERTY);
        if (filter == null || filter.isEmpty()) {
            return true;
        }
        final String parameterName = parameter == null ? null : parameter.getParameterName();
        return parameterName != null && parameterName.contains(filter);
    }

    private void emitPerBranchGradientDebug(double[] analyticGradient,
                                            double[][] perNodeGradient,
                                            double[][] perNodeNumericalGradient,
                                            int[] destinationIndexByNode,
                                            int dimGradient) {
        final double[] savedValues = parameter.getParameterValues();
        final double[] numericGradient = NumericalDerivative.gradient(numeric1, savedValues);
        for (int i = 0; i < savedValues.length; ++i) {
            parameter.setParameterValue(i, savedValues[i]);
        }
        treeDataLikelihood.makeDirty();

        StringBuilder sb = new StringBuilder();
        sb.append("Branch-specific gradient diagnostics (")
                .append(BRANCH_DEBUG_PROPERTY)
                .append("=true")
                .append(")\n");

        final double[] localNumericSum = new double[analyticGradient.length];
        boolean hasLocalNumeric = false;

        if (analyticGradient.length > 0) {
            final int[] sweepIndices;
            if (analyticGradient.length <= 8) {
                sweepIndices = new int[analyticGradient.length];
                for (int i = 0; i < analyticGradient.length; ++i) {
                    sweepIndices[i] = i;
                }
            } else if (analyticGradient.length > 9) {
                sweepIndices = new int[]{0, 6, 7, 8, 9};
            } else {
                sweepIndices = new int[]{0};
            }
            final double[] steps = new double[]{1e-6, 1e-5, 1e-4, 1e-3};
            for (final int index : sweepIndices) {
                if (index < 0 || index >= analyticGradient.length) {
                    continue;
                }
                sb.append("index").append(index).append(" finite-difference sweep:");
                for (double h : steps) {
                    final double slope = finiteDifferenceAtIndex(savedValues, index, h);
                    sb.append(" h=").append(h).append(" slope=").append(slope);
                }
                sb.append('\n');
            }
        }

        for (int nodeIndex = 0; nodeIndex < tree.getNodeCount(); ++nodeIndex) {
            double[] nodeGradient = perNodeGradient[nodeIndex];
            if (nodeGradient == null) {
                continue;
            }
            final double[] nodeNumericGradient = perNodeNumericalGradient[nodeIndex];
            NodeRef node = tree.getNode(nodeIndex);
            int destinationIndex = destinationIndexByNode[nodeIndex];

            double[] analyticSlice = new double[dimGradient];
            double[] numericSlice = new double[dimGradient];
            for (int j = 0; j < dimGradient; ++j) {
                final int gradientIndex = getResultIndex(destinationIndex, j, dimGradient, analyticGradient.length);
                analyticSlice[j] = analyticGradient[gradientIndex];
                numericSlice[j] = numericGradient[gradientIndex];
                if (nodeNumericGradient != null
                        && gradientIndex < localNumericSum.length
                        && j < nodeNumericGradient.length) {
                    localNumericSum[gradientIndex] += nodeNumericGradient[j];
                    hasLocalNumeric = true;
                }
            }

            sb.append("node=")
                    .append(node.getNumber())
                    .append(" isRoot=")
                    .append(tree.isRoot(node))
                    .append(" paramIndex=")
                    .append(destinationIndex)
                    .append(" localAnalytic=")
                    .append(new dr.math.matrixAlgebra.Vector(nodeGradient))
                    .append(" localNumeric=")
                    .append(nodeNumericGradient == null
                            ? "<unavailable>"
                            : new dr.math.matrixAlgebra.Vector(nodeNumericGradient))
                    .append(" analyticSlice=")
                    .append(new dr.math.matrixAlgebra.Vector(analyticSlice))
                    .append(" numericSlice=")
                    .append(new dr.math.matrixAlgebra.Vector(numericSlice))
                    .append('\n');
        }

        if (hasLocalNumeric) {
            sb.append("localNumericSum=")
                    .append(new dr.math.matrixAlgebra.Vector(localNumericSum))
                    .append('\n');
            sb.append("analyticMinusLocalNumericSum=")
                    .append(new dr.math.matrixAlgebra.Vector(vectorDifference(analyticGradient, localNumericSum)))
                    .append('\n');
            sb.append("globalNumericMinusLocalNumericSum=")
                    .append(new dr.math.matrixAlgebra.Vector(vectorDifference(numericGradient, localNumericSum)))
                    .append('\n');
        }

        System.err.print(sb);
    }

    private static double[] vectorDifference(final double[] left, final double[] right) {
        final int length = Math.min(left.length, right.length);
        final double[] out = new double[length];
        for (int i = 0; i < length; ++i) {
            out[i] = left[i] - right[i];
        }
        return out;
    }

    private double finiteDifferenceAtIndex(double[] baseValues, int index, double h) {
        final double old = baseValues[index];

        baseValues[index] = old + h;
        parameter.setParameterValue(index, baseValues[index]);
        treeDataLikelihood.makeDirty();
        final double plus = treeDataLikelihood.getLogLikelihood();

        baseValues[index] = old - h;
        parameter.setParameterValue(index, baseValues[index]);
        treeDataLikelihood.makeDirty();
        final double minus = treeDataLikelihood.getLogLikelihood();

        baseValues[index] = old;
        parameter.setParameterValue(index, old);
        treeDataLikelihood.makeDirty();

        return (plus - minus) / (2.0 * h);
    }

    private int getParameterIndexFromNode(NodeRef node) {
        return branchProvider.getParameterIndexFromNode(node);
//        return (branchRateModel == null) ? node.getNumber() : branchRateModel.getParameterIndexFromNode(node);
    }

    public List<ContinuousTraitGradientForBranch.ContinuousProcessParameterGradient.DerivationParameter> getDerivationParameter() {
        return ((ContinuousTraitGradientForBranch.ContinuousProcessParameterGradient) branchProvider).getDerivationParameter();
    }

    public List<ContinuousTraitGradientForBranch.ContinuousProcessParameterGradient.DerivationParameter> getDerivationParameterOrNull() {
        if (branchProvider instanceof ContinuousTraitGradientForBranch.ContinuousProcessParameterGradient) {
            return ((ContinuousTraitGradientForBranch.ContinuousProcessParameterGradient) branchProvider).getDerivationParameter();
        }
        return null;
    }

    public List<ContinuousTraitGradientForBranch.ProcessGradientTarget> getProcessGradientTargetsOrNull() {
        if (branchProvider instanceof ContinuousTraitGradientForBranch.TargetBasedContinuousProcessParameterGradient) {
            return ((ContinuousTraitGradientForBranch.TargetBasedContinuousProcessParameterGradient) branchProvider).getTargets();
        }
        return null;
    }

    private MultivariateFunction numeric1 = new MultivariateFunction() {
        @Override
        public double evaluate(double[] argument) {

            for (int i = 0; i < argument.length; ++i) {
                parameter.setParameterValue(i, argument[i]);
            }

            treeDataLikelihood.makeDirty();
            return treeDataLikelihood.getLogLikelihood();
        }

        @Override
        public int getNumArguments() {
            return parameter.getDimension();
        }

        @Override
        public double getLowerBound(int n) {
            return 0;
        }

        @Override
        public double getUpperBound(int n) {
            return Double.POSITIVE_INFINITY;
        }
    };

    /**
     * Debug helper: finite-difference gradient of this branch-specific provider at current parameter state.
     * Restores parameter values and dirties the likelihood before returning.
     */
    public double[] getNumericalGradientLogDensityDebug() {
        final double[] savedValues = parameter.getParameterValues();
        final double[] numeric = NumericalDerivative.gradient(numeric1, savedValues);
        for (int i = 0; i < savedValues.length; ++i) {
            parameter.setParameterValue(i, savedValues[i]);
        }
        treeDataLikelihood.makeDirty();
        return numeric;
    }

    /**
     * Snapshot of branch/global gradient consistency at the current state.
     * <p>
     * This is designed for regression tests and replay debugging so we can
     * detect stale/global-local mismatches without rebuilding ad-hoc tooling.
     */
    public static final class GradientConsistencySnapshot {
        public final double[] fullAnalyticFirst;
        public final double[] fullAnalyticSecond;
        public final double[] fullNumeric;
        public final double[] localAnalyticSum;
        public final double[] localNumericSum;

        private GradientConsistencySnapshot(double[] fullAnalyticFirst,
                                            double[] fullAnalyticSecond,
                                            double[] fullNumeric,
                                            double[] localAnalyticSum,
                                            double[] localNumericSum) {
            this.fullAnalyticFirst = fullAnalyticFirst;
            this.fullAnalyticSecond = fullAnalyticSecond;
            this.fullNumeric = fullNumeric;
            this.localAnalyticSum = localAnalyticSum;
            this.localNumericSum = localNumericSum;
        }
    }

    public GradientConsistencySnapshot computeGradientConsistencySnapshot() {
        final int dimGradient = branchProvider.getDimension();
        final int parameterDimension = parameter.getDimension();

        final double[] fullAnalyticFirst = getGradientLogDensity();
        final double[] fullNumeric = getNumericalGradientLogDensityDebug();
        final double[] localAnalyticSum = new double[parameterDimension];
        final double[] localNumericSum = new double[parameterDimension];

        for (int i = 0; i < tree.getNodeCount(); ++i) {
            final NodeRef node = tree.getNode(i);
            if (tree.isRoot(node)) {
                continue;
            }
            final List<BranchSufficientStatistics> statisticsForNode = treeTraitProvider.getTrait(tree, node);
            final BranchSufficientStatistics statistics = statisticsForNode.get(0);
            final double[] localAnalytic = branchProvider.getGradientForBranch(statistics, node);
            final double[] localNumeric = branchProvider.getNumericalGradientForBranch(statistics, node);
            final int destinationIndex = getParameterIndexFromNode(node);

            for (int j = 0; j < dimGradient; ++j) {
                final int index = getResultIndex(destinationIndex, j, dimGradient, parameterDimension);
                localAnalyticSum[index] += localAnalytic[j];
                if (localNumeric != null && j < localNumeric.length) {
                    localNumericSum[index] += localNumeric[j];
                }
            }
        }

        final double[] fullAnalyticSecond = getGradientLogDensity();
        return new GradientConsistencySnapshot(
                fullAnalyticFirst,
                fullAnalyticSecond,
                fullNumeric,
                localAnalyticSum,
                localNumericSum
        );
    }

    /**
     * Debug helper: compare two immediate consecutive branch-stat passes at fixed state.
     * Intended for replay diagnosis of stale/cached branch sufficient statistics.
     */
    public String getConsecutivePassStatisticsDiffReport() {
        final List<NodeStatsDigest> first = collectNodeStatsDigest();
        final List<NodeStatsDigest> second = collectNodeStatsDigest();
        final StringBuilder sb = new StringBuilder(4096);
        sb.append("consecutivePassStatsDiff\n");
        double globalMax = 0.0;
        int globalMaxNode = -1;
        String globalMaxField = "<none>";
        int changedNodes = 0;
        for (int i = 0; i < first.size(); ++i) {
            final NodeStatsDigest a = first.get(i);
            final NodeStatsDigest b = second.get(i);
            if (a.nodeNumber != b.nodeNumber) {
                throw new IllegalStateException("Node ordering mismatch across consecutive digest passes.");
            }
            final double dBranchActualization = maxAbsDiff(a.branchActualization, b.branchActualization);
            final double dBranchDisplacement = maxAbsDiff(a.branchDisplacement, b.branchDisplacement);
            final double dBranchPrecision = maxAbsDiff(a.branchPrecision, b.branchPrecision);
            final double dAboveMean = maxAbsDiff(a.aboveMean, b.aboveMean);
            final double dAbovePrecision = maxAbsDiff(a.abovePrecision, b.abovePrecision);
            final double dBelowMean = maxAbsDiff(a.belowMean, b.belowMean);
            final double dBelowPrecision = maxAbsDiff(a.belowPrecision, b.belowPrecision);
            final double dLocalAnalytic = maxAbsDiff(a.localAnalytic, b.localAnalytic);

            final double nodeMax = max(
                    dBranchActualization, dBranchDisplacement, dBranchPrecision,
                    dAboveMean, dAbovePrecision, dBelowMean, dBelowPrecision, dLocalAnalytic
            );
            if (nodeMax > 0.0) {
                changedNodes++;
            }
            if (nodeMax > globalMax) {
                globalMax = nodeMax;
                globalMaxNode = a.nodeNumber;
                globalMaxField = argmaxField(
                        dBranchActualization, dBranchDisplacement, dBranchPrecision,
                        dAboveMean, dAbovePrecision, dBelowMean, dBelowPrecision, dLocalAnalytic
                );
            }
            if (nodeMax > 0.0) {
                sb.append("\tnode=").append(a.nodeNumber)
                        .append(" isExternal=").append(a.isExternal)
                        .append(" destinationIndex=").append(a.destinationIndex)
                        .append(" max=").append(nodeMax)
                        .append(" dBranchActualization=").append(dBranchActualization)
                        .append(" dBranchDisplacement=").append(dBranchDisplacement)
                        .append(" dBranchPrecision=").append(dBranchPrecision)
                        .append(" dAboveMean=").append(dAboveMean)
                        .append(" dAbovePrecision=").append(dAbovePrecision)
                        .append(" dBelowMean=").append(dBelowMean)
                        .append(" dBelowPrecision=").append(dBelowPrecision)
                        .append(" dLocalAnalytic=").append(dLocalAnalytic)
                        .append('\n');
            }
        }
        sb.append("\tchangedNodes=").append(changedNodes)
                .append(" totalNodes=").append(first.size())
                .append(" globalMax=").append(globalMax)
                .append(" globalMaxNode=").append(globalMaxNode)
                .append(" globalMaxField=").append(globalMaxField)
                .append('\n');
        return sb.toString();
    }

    private List<NodeStatsDigest> collectNodeStatsDigest() {
        final List<NodeStatsDigest> out = new ArrayList<NodeStatsDigest>();
        for (int i = 0; i < tree.getNodeCount(); ++i) {
            final NodeRef node = tree.getNode(i);
            if (tree.isRoot(node)) {
                continue;
            }
            final List<BranchSufficientStatistics> statisticsForNode = treeTraitProvider.getTrait(tree, node);
            final BranchSufficientStatistics statistics = statisticsForNode.get(0);
            final int destinationIndex = getParameterIndexFromNode(node);
            final double[] localAnalytic = branchProvider.getGradientForBranch(statistics, node);
            out.add(new NodeStatsDigest(
                    node.getNumber(),
                    tree.isExternal(node),
                    destinationIndex,
                    copyDenseMatrix(statistics.getBranch().getRawActualization()),
                    copyDenseVectorColumn(statistics.getBranch().getRawDisplacement()),
                    copyDenseMatrix(statistics.getBranch().getRawPrecision()),
                    copyDenseVectorColumn(statistics.getAbove().getRawMean()),
                    copyDenseMatrix(statistics.getAbove().getRawPrecision()),
                    copyDenseVectorColumn(statistics.getBelow().getRawMean()),
                    copyDenseMatrix(statistics.getBelow().getRawPrecision()),
                    Arrays.copyOf(localAnalytic, localAnalytic.length)
            ));
        }
        return out;
    }

    private static double[] copyDenseMatrix(final DenseMatrix64F matrix) {
        final int length = matrix.getNumRows() * matrix.getNumCols();
        final double[] out = new double[length];
        int k = 0;
        for (int i = 0; i < matrix.getNumRows(); ++i) {
            for (int j = 0; j < matrix.getNumCols(); ++j) {
                out[k++] = matrix.get(i, j);
            }
        }
        return out;
    }

    private static double[] copyDenseVectorColumn(final DenseMatrix64F vector) {
        final double[] out = new double[vector.getNumRows()];
        for (int i = 0; i < vector.getNumRows(); ++i) {
            out[i] = vector.get(i, 0);
        }
        return out;
    }

    private static double maxAbsDiff(final double[] left, final double[] right) {
        final int length = Math.min(left.length, right.length);
        double max = 0.0;
        for (int i = 0; i < length; ++i) {
            max = Math.max(max, Math.abs(left[i] - right[i]));
        }
        return max;
    }

    private static double max(final double... values) {
        double m = Double.NEGATIVE_INFINITY;
        for (double v : values) {
            m = Math.max(m, v);
        }
        return m;
    }

    private static String argmaxField(final double dBranchActualization,
                                      final double dBranchDisplacement,
                                      final double dBranchPrecision,
                                      final double dAboveMean,
                                      final double dAbovePrecision,
                                      final double dBelowMean,
                                      final double dBelowPrecision,
                                      final double dLocalAnalytic) {
        String best = "dBranchActualization";
        double max = dBranchActualization;
        if (dBranchDisplacement > max) { max = dBranchDisplacement; best = "dBranchDisplacement"; }
        if (dBranchPrecision > max) { max = dBranchPrecision; best = "dBranchPrecision"; }
        if (dAboveMean > max) { max = dAboveMean; best = "dAboveMean"; }
        if (dAbovePrecision > max) { max = dAbovePrecision; best = "dAbovePrecision"; }
        if (dBelowMean > max) { max = dBelowMean; best = "dBelowMean"; }
        if (dBelowPrecision > max) { max = dBelowPrecision; best = "dBelowPrecision"; }
        if (dLocalAnalytic > max) { best = "dLocalAnalytic"; }
        return best;
    }

    private static final class NodeStatsDigest {
        private final int nodeNumber;
        private final boolean isExternal;
        private final int destinationIndex;
        private final double[] branchActualization;
        private final double[] branchDisplacement;
        private final double[] branchPrecision;
        private final double[] aboveMean;
        private final double[] abovePrecision;
        private final double[] belowMean;
        private final double[] belowPrecision;
        private final double[] localAnalytic;

        private NodeStatsDigest(final int nodeNumber,
                                final boolean isExternal,
                                final int destinationIndex,
                                final double[] branchActualization,
                                final double[] branchDisplacement,
                                final double[] branchPrecision,
                                final double[] aboveMean,
                                final double[] abovePrecision,
                                final double[] belowMean,
                                final double[] belowPrecision,
                                final double[] localAnalytic) {
            this.nodeNumber = nodeNumber;
            this.isExternal = isExternal;
            this.destinationIndex = destinationIndex;
            this.branchActualization = branchActualization;
            this.branchDisplacement = branchDisplacement;
            this.branchPrecision = branchPrecision;
            this.aboveMean = aboveMean;
            this.abovePrecision = abovePrecision;
            this.belowMean = belowMean;
            this.belowPrecision = belowPrecision;
            this.localAnalytic = localAnalytic;
        }
    }

    /**
     * One-shot JSON dump at current state for Python side-by-side replay debugging.
     */
    public void dumpCurrentBranchStateForPythonDebug(final String path,
                                                     final double[] fullAnalyticGradient,
                                                     final double[] fullNumericGradient) {
        final boolean includeLocalNumeric =
                Boolean.getBoolean("beast.debug.branchDump.includeLocalNumeric");
        final StringBuilder sb = new StringBuilder(64 * 1024);
        sb.append("{\n");
        sb.append("\"parameterName\":\"").append(escapeJson(parameter.getParameterName())).append("\",\n");
        sb.append("\"parameterDimension\":").append(parameter.getDimension()).append(",\n");
        sb.append("\"parameterValues\":").append(jsonArray(parameter.getParameterValues())).append(",\n");
        sb.append("\"localDimension\":").append(branchProvider.getDimension()).append(",\n");
        sb.append("\"fullAnalytic\":").append(jsonArray(fullAnalyticGradient)).append(",\n");
        sb.append("\"fullNumeric\":").append(jsonArray(fullNumericGradient)).append(",\n");
        sb.append("\"branches\":[\n");

        boolean first = true;
        for (int i = 0; i < tree.getNodeCount(); ++i) {
            final NodeRef node = tree.getNode(i);
            if (tree.isRoot(node)) {
                continue;
            }
            final List<BranchSufficientStatistics> statisticsForNode = treeTraitProvider.getTrait(tree, node);
            final BranchSufficientStatistics statistics = statisticsForNode.get(0);
            final double[] localAnalytic = branchProvider.getGradientForBranch(statistics, node);
            final double[] localNumeric = includeLocalNumeric
                    ? branchProvider.getNumericalGradientForBranch(statistics, node)
                    : null;
            final NormalSufficientStatistics below = statistics.getBelow();
            final NormalSufficientStatistics above = statistics.getAbove();
            final NormalSufficientStatistics jointStatistics =
                    BranchRateGradient.ContinuousTraitGradientForBranch.Dense.computeJointStatistics(
                            below, above, treeDataLikelihood.getDataLikelihoodDelegate().getTraitDim());
            final int dim = treeDataLikelihood.getDataLikelihoodDelegate().getTraitDim();
            final DenseMatrix64F delta = new DenseMatrix64F(dim, 1);
            for (int row = 0; row < dim; ++row) {
                delta.unsafe_set(row, 0, jointStatistics.getMean(row) - above.getMean(row));
            }
            final DenseMatrix64F gradN = new DenseMatrix64F(dim, 1);
            CommonOps.multTransA(above.getRawPrecision(), delta, gradN);
            final NodeRef parent = tree.getParent(node);
            final boolean parentIsRoot = parent != null && tree.isRoot(parent);

            if (!first) {
                sb.append(",\n");
            }
            first = false;
            sb.append("{");
            sb.append("\"node\":").append(node.getNumber()).append(",");
            sb.append("\"isExternal\":").append(tree.isExternal(node)).append(",");
            sb.append("\"parentIsRoot\":").append(parentIsRoot).append(",");
            sb.append("\"destinationIndex\":").append(getParameterIndexFromNode(node)).append(",");
            sb.append("\"branchActualization\":").append(jsonMatrix(statistics.getBranch().getRawActualization())).append(",");
            sb.append("\"branchDisplacement\":").append(jsonVectorFromDense(statistics.getBranch().getRawDisplacement())).append(",");
            sb.append("\"branchPrecision\":").append(jsonMatrix(statistics.getBranch().getRawPrecision())).append(",");
            sb.append("\"aboveMean\":").append(jsonVectorFromDense(statistics.getAbove().getRawMean())).append(",");
            sb.append("\"abovePrecision\":").append(jsonMatrix(statistics.getAbove().getRawPrecision())).append(",");
            sb.append("\"belowMean\":").append(jsonVectorFromDense(statistics.getBelow().getRawMean())).append(",");
            sb.append("\"belowPrecision\":").append(jsonMatrix(statistics.getBelow().getRawPrecision())).append(",");
            sb.append("\"gradN\":").append(jsonVectorFromDense(gradN)).append(",");
            sb.append("\"missingIndices\":").append(jsonIntArray(statistics.getMissing())).append(",");
            sb.append("\"localAnalytic\":").append(jsonArray(localAnalytic)).append(",");
            sb.append("\"localNumeric\":").append(localNumeric == null ? "null" : jsonArray(localNumeric));
            sb.append("}");
        }

        final NodeRef root = tree.getRoot();
        final List<BranchSufficientStatistics> rootStatisticsForNode = treeTraitProvider.getTrait(tree, root);
        final BranchSufficientStatistics rootStatistics = rootStatisticsForNode.get(0);
        final double[] rootLocalAnalytic = branchProvider.getGradientForBranch(rootStatistics, root);
        final double[] rootLocalNumeric = includeLocalNumeric
                ? branchProvider.getNumericalGradientForBranch(rootStatistics, root)
                : null;
        if (!first) {
            sb.append(",\n");
        }
        sb.append("{");
        sb.append("\"node\":").append(root.getNumber()).append(",");
        sb.append("\"isRoot\":true,");
        sb.append("\"isExternal\":").append(tree.isExternal(root)).append(",");
        sb.append("\"destinationIndex\":").append(getParameterIndexFromNode(root)).append(",");
        sb.append("\"aboveMean\":").append(jsonVectorFromDense(rootStatistics.getAbove().getRawMean())).append(",");
        sb.append("\"abovePrecision\":").append(jsonMatrix(rootStatistics.getAbove().getRawPrecision())).append(",");
        sb.append("\"belowMean\":").append(jsonVectorFromDense(rootStatistics.getBelow().getRawMean())).append(",");
        sb.append("\"belowPrecision\":").append(jsonMatrix(rootStatistics.getBelow().getRawPrecision())).append(",");
        sb.append("\"missingIndices\":").append(jsonIntArray(rootStatistics.getMissing())).append(",");
        sb.append("\"localAnalytic\":").append(jsonArray(rootLocalAnalytic)).append(",");
        sb.append("\"localNumeric\":").append(rootLocalNumeric == null ? "null" : jsonArray(rootLocalNumeric));
        sb.append("}");
        sb.append("\n]\n}\n");

        try {
            Files.write(Paths.get(path), sb.toString().getBytes(StandardCharsets.UTF_8));
        } catch (Exception e) {
            throw new RuntimeException("Failed writing branch-state debug dump to " + path, e);
        }
    }

    private static String jsonArray(final double[] values) {
        final StringBuilder sb = new StringBuilder();
        sb.append("[");
        for (int i = 0; i < values.length; ++i) {
            if (i > 0) {
                sb.append(",");
            }
            sb.append(values[i]);
        }
        sb.append("]");
        return sb.toString();
    }

    private static String jsonIntArray(final int[] values) {
        final StringBuilder sb = new StringBuilder();
        sb.append("[");
        for (int i = 0; i < values.length; ++i) {
            if (i > 0) {
                sb.append(",");
            }
            sb.append(values[i]);
        }
        sb.append("]");
        return sb.toString();
    }

    private static String jsonVectorFromDense(final DenseMatrix64F values) {
        final StringBuilder sb = new StringBuilder();
        sb.append("[");
        for (int i = 0; i < values.getNumRows(); ++i) {
            if (i > 0) {
                sb.append(",");
            }
            sb.append(values.get(i, 0));
        }
        sb.append("]");
        return sb.toString();
    }

    private static String jsonMatrix(final DenseMatrix64F values) {
        final StringBuilder sb = new StringBuilder();
        sb.append("[");
        for (int i = 0; i < values.getNumRows(); ++i) {
            if (i > 0) {
                sb.append(",");
            }
            sb.append("[");
            for (int j = 0; j < values.getNumCols(); ++j) {
                if (j > 0) {
                    sb.append(",");
                }
                sb.append(values.get(i, j));
            }
            sb.append("]");
        }
        sb.append("]");
        return sb.toString();
    }

    private static String escapeJson(final String value) {
        if (value == null) {
            return "null";
        }
        return value.replace("\\", "\\\\").replace("\"", "\\\"");
    }

    @Override
    public String getReport() {

        double[] savedValues = parameter.getParameterValues();
        double[] testGradient = NumericalDerivative.gradient(numeric1, parameter.getParameterValues());
        for (int i = 0; i < savedValues.length; ++i) {
            parameter.setParameterValue(i, savedValues[i]);
        }

        StringBuilder sb = new StringBuilder();
        sb.append("Peeling: ").append(new dr.math.matrixAlgebra.Vector(getGradientLogDensity()));
        sb.append("\n");
        sb.append("numeric: ").append(new dr.math.matrixAlgebra.Vector(testGradient));
        sb.append("\n");

        return sb.toString();
    }

    private static final boolean DEBUG = false;

    @Override
    public LogColumn[] getColumns() {

        LogColumn[] columns = new LogColumn[1];
        columns[0] = new LogColumn.Default("gradient report", new Object() {
            @Override
            public String toString() {
                return "\n" + getReport();
            }
        });

        return columns;
    }
}
