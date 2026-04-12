/*
 * DiffusionParametersGradient.java
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

package dr.evomodel.treedatalikelihood.hmc;

import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.evomodel.treedatalikelihood.continuous.BranchSpecificGradient;
import dr.evomodel.treedatalikelihood.continuous.ContinuousTraitGradientForBranch;
import dr.inference.hmc.CompoundGradient;
import dr.inference.hmc.GradientWrtParameterProvider;
import dr.inference.model.CompoundParameter;
import dr.inference.model.Likelihood;
import dr.inference.model.Parameter;
import dr.xml.Reportable;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.IdentityHashMap;
import java.util.List;
import java.util.Map;
import java.util.concurrent.atomic.AtomicLong;

public class DiffusionParametersGradient implements GradientWrtParameterProvider, Reportable {
    private static final AtomicLong DEBUG_CALL_COUNTER = new AtomicLong(0L);

    private final TreeDataLikelihood likelihood;
    private final int dim;
    private final BranchSpecificGradient branchSpecificGradient;
    private final CompoundParameter compoundParameter;
    private final CompoundGradient parametersGradients;
    private final List<GradientSlice> gradientSlices;

    public DiffusionParametersGradient(BranchSpecificGradient branchSpecificGradient, CompoundGradient parametersGradients) {

        this.branchSpecificGradient = branchSpecificGradient;
        this.likelihood = (TreeDataLikelihood) branchSpecificGradient.getLikelihood();

        compoundParameter = new CompoundParameter(null);
        this.parametersGradients = parametersGradients;
        this.gradientSlices = buildGradientSlices(parametersGradients, compoundParameter);
        this.dim = sumResultDimensions(gradientSlices);
        validateDerivationOrder();

        if (Boolean.getBoolean("beast.debug.compoundParameter.identitySummary")) {
            emitCompoundParameterIdentitySummary("DiffusionParametersGradient", compoundParameter);
        }

    }

    private List<GradientSlice> buildGradientSlices(CompoundGradient parametersGradients, CompoundParameter parameter) {
        final List<GradientSlice> slices = new ArrayList<GradientSlice>();
        int sourceOffset = 0;
        int resultOffset = 0;
        int dimTrait = likelihood.getDataLikelihoodDelegate().getTraitDim();
        for (GradientWrtParameterProvider gradient : parametersGradients.getDerivativeList()) {
            assert gradient instanceof AbstractDiffusionGradient : "Gradients must all be instances of AbstractDiffusionGradient.";
            final AbstractDiffusionGradient diffusionGradient = (AbstractDiffusionGradient) gradient;
            parameter.addParameter(gradient.getParameter());
            final int sourceDimension = diffusionGradient.getDerivationParameter().getDimension(dimTrait);
            final int resultDimension = gradient.getDimension();
            slices.add(new GradientSlice(diffusionGradient, sourceOffset, sourceDimension, resultOffset, resultDimension));
            sourceOffset += sourceDimension;
            resultOffset += resultDimension;
        }
        return Collections.unmodifiableList(slices);
    }

    private int sumResultDimensions(List<GradientSlice> slices) {
        int total = 0;
        for (GradientSlice slice : slices) {
            total += slice.resultDimension;
        }
        return total;
    }

    private void validateDerivationOrder() {
        final List<ContinuousTraitGradientForBranch.ContinuousProcessParameterGradient.DerivationParameter> branchDerivationOrder =
                branchSpecificGradient.getDerivationParameterOrNull();
        if (branchDerivationOrder != null) {
            validateAgainstDerivationOrder(branchDerivationOrder);
            return;
        }

        final List<ContinuousTraitGradientForBranch.ProcessGradientTarget> branchTargetOrder =
                branchSpecificGradient.getProcessGradientTargetsOrNull();
        if (branchTargetOrder != null) {
            validateAgainstTargetOrder(branchTargetOrder);
            return;
        }

        throw new IllegalArgumentException(
                "Unsupported branch gradient provider for diffusion-parameter validation: unable to discover derivation or target order.");
    }

    private void validateAgainstDerivationOrder(
            final List<ContinuousTraitGradientForBranch.ContinuousProcessParameterGradient.DerivationParameter> branchOrder) {
        if (branchOrder.size() != gradientSlices.size()) {
            throw new IllegalArgumentException("Diffusion gradient block count does not match branch derivation count.");
        }
        final int dimTrait = likelihood.getDataLikelihoodDelegate().getTraitDim();
        for (int i = 0; i < gradientSlices.size(); i++) {
            final GradientSlice slice = gradientSlices.get(i);
            final ContinuousTraitGradientForBranch.ContinuousProcessParameterGradient.DerivationParameter expected = branchOrder.get(i);
            final ContinuousTraitGradientForBranch.ContinuousProcessParameterGradient.DerivationParameter actual =
                    slice.gradientProvider.getDerivationParameter();
            if (expected != actual) {
                throw new IllegalArgumentException("Diffusion gradient order mismatch: expected " + expected + " but found " + actual + ".");
            }
            final int expectedSourceDimension = expected.getDimension(dimTrait);
            if (slice.sourceDimension != expectedSourceDimension) {
                throw new IllegalArgumentException(
                        "Diffusion source dimension mismatch for " + actual + ": expected " + expectedSourceDimension +
                                " but found " + slice.sourceDimension + "."
                );
            }
        }
    }

    private void validateAgainstTargetOrder(final List<ContinuousTraitGradientForBranch.ProcessGradientTarget> branchOrder) {
        if (branchOrder.size() != gradientSlices.size()) {
            throw new IllegalArgumentException("Diffusion gradient block count does not match branch target count.");
        }
        final int dimTrait = likelihood.getDataLikelihoodDelegate().getTraitDim();
        for (int i = 0; i < gradientSlices.size(); i++) {
            final GradientSlice slice = gradientSlices.get(i);
            final ContinuousTraitGradientForBranch.ProcessGradientTarget expected = branchOrder.get(i);
            final ContinuousTraitGradientForBranch.ProcessGradientTarget actual = slice.gradientProvider.getProcessGradientTarget();
            if (!expected.getClass().equals(actual.getClass())) {
                throw new IllegalArgumentException(
                        "Diffusion gradient target order mismatch at index " + i +
                                ": expected " + expected.getClass().getSimpleName() +
                                " but found " + actual.getClass().getSimpleName() + "."
                );
            }
            final int expectedSourceDimension = expected.getDimension(dimTrait);
            final int actualSourceDimension = actual.getDimension(dimTrait);
            if (slice.sourceDimension != expectedSourceDimension || expectedSourceDimension != actualSourceDimension) {
                throw new IllegalArgumentException(
                        "Diffusion source dimension mismatch for target order at index " + i +
                                ": expected=" + expectedSourceDimension +
                                ", actual=" + actualSourceDimension +
                                ", slice=" + slice.sourceDimension + "."
                );
            }
        }
    }

    @Override
    public Likelihood getLikelihood() {
        return likelihood;
    }

    @Override
    public Parameter getParameter() {
        return compoundParameter;
    }

    @Override
    public int getDimension() {
        return dim;
    }

    @Override
    public double[] getGradientLogDensity() {
        if (Boolean.getBoolean("beast.gradientHook.diffusionParametersGradient.forceDirtyBeforeAnalytic")) {
            likelihood.makeDirty();
            if (Boolean.getBoolean("beast.gradientHook.diffusionParametersGradient.forceLikelihoodEval")) {
                likelihood.getLogLikelihood();
            }
        }
        double[] gradient = branchSpecificGradient.getGradientLogDensity();
        maybeEmitRecomputeConsistencyDebug(gradient);
        return getGradientLogDensity(gradient);
    }

    private double[] getGradientLogDensity(double[] gradient) {
        double[] result = new double[dim];
        for (GradientSlice slice : gradientSlices) {
            System.arraycopy(
                    slice.gradientProvider.getGradientLogDensity(gradient, slice.sourceOffset),
                    0,
                    result,
                    slice.resultOffset,
                    slice.resultDimension
            );
        }
        maybeEmitSliceDebug(gradient, result);
        return result;
    }

    private void maybeEmitSliceDebug(final double[] fullSourceGradient,
                                     final double[] slicedResultGradient) {
        if (!Boolean.getBoolean("beast.debug.diffusionSlice")) {
            return;
        }
        final long call = DEBUG_CALL_COUNTER.incrementAndGet();
        final int stride = Integer.getInteger("beast.debug.diffusionSliceStride", 1);
        if (stride <= 0 || (call % stride != 0)) {
            return;
        }

        final int[] focus = parseFocusIndices();
        final String focusParamContains = System.getProperty(
                "beast.debug.diffusionSlice.onlyParameterContains",
                "OrthogonalBlockDiagonalPolarStableMatrixParameter.native"
        );
        final double[] numericFull = branchSpecificGradient.getNumericalGradientLogDensityDebug();

        final StringBuilder sb = new StringBuilder();
        sb.append("diffusionSliceDebug call=").append(call).append('\n');
        sb.append("\tfullSource(analytic)");
        appendIndexValues(sb, fullSourceGradient, focus);
        sb.append('\n');
        sb.append("\tfullSource(numeric )");
        appendIndexValues(sb, numericFull, focus);
        sb.append('\n');
        sb.append("\tfullSource(absDiff)");
        appendIndexDiffs(sb, fullSourceGradient, numericFull, focus);
        sb.append('\n');

        for (GradientSlice slice : gradientSlices) {
            final String rawName = slice.gradientProvider.getRawParameter() == null
                    ? "<null>"
                    : slice.gradientProvider.getRawParameter().getParameterName();
            if (focusParamContains != null && !focusParamContains.isEmpty() && !rawName.contains(focusParamContains)) {
                continue;
            }
            final double[] analyticSlice = slice.gradientProvider.getGradientLogDensity(fullSourceGradient, slice.sourceOffset);
            final double[] numericSlice = slice.gradientProvider.getGradientLogDensity(numericFull, slice.sourceOffset);
            sb.append("\tslice parameter=").append(rawName)
                    .append(" sourceOffset=").append(slice.sourceOffset)
                    .append(" sourceDim=").append(slice.sourceDimension)
                    .append(" resultOffset=").append(slice.resultOffset)
                    .append(" resultDim=").append(slice.resultDimension)
                    .append('\n');
            sb.append("\t\tanalyticSlice=").append(Arrays.toString(analyticSlice)).append('\n');
            sb.append("\t\tnumericSlice =").append(Arrays.toString(numericSlice)).append('\n');
            sb.append("\t\tabsDiffSlice =").append(Arrays.toString(absDiff(analyticSlice, numericSlice))).append('\n');
        }

        System.err.print(sb.toString());
    }

    private int[] parseFocusIndices() {
        final String raw = System.getProperty("beast.debug.diffusionSlice.indices", "6,7,8");
        final String[] tokens = raw.split(",");
        final List<Integer> parsed = new ArrayList<>();
        for (String token : tokens) {
            final String t = token.trim();
            if (t.isEmpty()) {
                continue;
            }
            try {
                parsed.add(Integer.parseInt(t));
            } catch (NumberFormatException ignored) {
                // debug-only parsing
            }
        }
        if (parsed.isEmpty()) {
            parsed.add(6);
            parsed.add(7);
            parsed.add(8);
        }
        final int[] out = new int[parsed.size()];
        for (int i = 0; i < parsed.size(); ++i) {
            out[i] = parsed.get(i);
        }
        return out;
    }

    private void appendIndexValues(final StringBuilder sb, final double[] vector, final int[] indices) {
        for (int idx : indices) {
            if (idx >= 0 && idx < vector.length) {
                sb.append(" i").append(idx).append('=').append(vector[idx]);
            } else {
                sb.append(" i").append(idx).append("=<oob>");
            }
        }
    }

    private void appendIndexDiffs(final StringBuilder sb, final double[] left, final double[] right, final int[] indices) {
        for (int idx : indices) {
            if (idx >= 0 && idx < left.length && idx < right.length) {
                sb.append(" i").append(idx).append('=').append(Math.abs(left[idx] - right[idx]));
            } else {
                sb.append(" i").append(idx).append("=<oob>");
            }
        }
    }

    private double[] absDiff(final double[] a, final double[] b) {
        final int n = Math.min(a.length, b.length);
        final double[] out = new double[n];
        for (int i = 0; i < n; ++i) {
            out[i] = Math.abs(a[i] - b[i]);
        }
        return out;
    }

    private void maybeEmitRecomputeConsistencyDebug(final double[] baselineGradient) {
        if (!Boolean.getBoolean("beast.debug.diffusionRecomputeConsistency")) {
            return;
        }
        final long call = DEBUG_CALL_COUNTER.incrementAndGet();
        final int stride = Integer.getInteger("beast.debug.diffusionRecomputeConsistencyStride", 1);
        if (stride <= 0 || (call % stride != 0)) {
            return;
        }

        likelihood.makeDirty();
        if (Boolean.getBoolean("beast.debug.diffusionRecomputeConsistencyForceLikelihood")) {
            likelihood.getLogLikelihood();
        }
        final double[] recomputedGradient = branchSpecificGradient.getGradientLogDensity();

        double maxAbs = 0.0;
        int maxIdx = -1;
        final int n = Math.min(baselineGradient.length, recomputedGradient.length);
        for (int i = 0; i < n; ++i) {
            final double diff = Math.abs(baselineGradient[i] - recomputedGradient[i]);
            if (diff > maxAbs) {
                maxAbs = diff;
                maxIdx = i;
            }
        }

        System.err.println("diffusionRecomputeConsistency call=" + call
                + " maxAbsDiff=" + maxAbs
                + " maxIdx=" + maxIdx
                + " baseline=" + (maxIdx >= 0 ? baselineGradient[maxIdx] : Double.NaN)
                + " recomputed=" + (maxIdx >= 0 ? recomputedGradient[maxIdx] : Double.NaN));
    }

    @Override
    public String getReport() {
        return "diffusionGradient." + compoundParameter.getParameterName() + "\n" +
                GradientWrtParameterProvider.getReportAndCheckForError(this,
                        Double.NEGATIVE_INFINITY, Double.POSITIVE_INFINITY,
                GradientWrtParameterProvider.TOLERANCE);
    }

    private static final class GradientSlice {
        private final AbstractDiffusionGradient gradientProvider;
        private final int sourceOffset;
        private final int sourceDimension;
        private final int resultOffset;
        private final int resultDimension;

        private GradientSlice(AbstractDiffusionGradient gradientProvider,
                              int sourceOffset,
                              int sourceDimension,
                              int resultOffset,
                              int resultDimension) {
            this.gradientProvider = gradientProvider;
            this.sourceOffset = sourceOffset;
            this.sourceDimension = sourceDimension;
            this.resultOffset = resultOffset;
            this.resultDimension = resultDimension;
        }
    }

    private static void emitCompoundParameterIdentitySummary(final String label,
                                                             final CompoundParameter compoundParameter) {
        final Map<Parameter, Integer> multiplicity = new IdentityHashMap<Parameter, Integer>();
        for (int i = 0; i < compoundParameter.getParameterCount(); ++i) {
            final Parameter child = compoundParameter.getParameter(i);
            final Integer previous = multiplicity.get(child);
            multiplicity.put(child, previous == null ? 1 : previous + 1);
        }
        int duplicatedChildRefs = 0;
        for (Integer count : multiplicity.values()) {
            if (count != null && count > 1) {
                duplicatedChildRefs++;
            }
        }
        System.err.println("compoundParameterIdentitySummary label=" + label
                + " id=" + System.identityHashCode(compoundParameter)
                + " parameterCount=" + compoundParameter.getParameterCount()
                + " identityUniqueCount=" + multiplicity.size()
                + " duplicatedChildRefs=" + duplicatedChildRefs
                + " dimension=" + compoundParameter.getDimension());
    }
}
