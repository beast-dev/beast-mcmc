/*
 * AbstractDiffusionGradient.java
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

import dr.evomodel.treedatalikelihood.continuous.BranchSpecificGradient;
import dr.evomodel.treedatalikelihood.continuous.ContinuousTraitGradientForBranch;
import dr.inference.model.AbstractBlockDiagonalTwoByTwoMatrixParameter;
import dr.inference.hmc.GradientWrtParameterProvider;
import dr.inference.model.DiagonalMatrix;
import dr.inference.model.Likelihood;
import dr.inference.model.MatrixParameterInterface;
import dr.inference.model.Parameter;
import dr.math.MultivariateFunction;
import dr.math.NumericalDerivative;
import dr.math.matrixAlgebra.Vector;
import dr.xml.Reportable;

import java.util.List;
import java.util.concurrent.atomic.AtomicLong;

public abstract class AbstractDiffusionGradient implements GradientWrtParameterProvider, Reportable {
    private static final double REPORT_SMALL_COMPONENT_THRESHOLD = 5.0e-3;
    private static final AtomicLong DEBUG_PARAMETER_DIFFUSION_CALL_COUNTER = new AtomicLong(0L);

    private final Likelihood likelihood;
    //    private final Parameter parameter;
    private final double lowerBound;
    private final double upperBound;
    protected int offset;

    AbstractDiffusionGradient(Likelihood likelihood, double upperBound, double lowerBound) {
        this.likelihood = likelihood;
//        this.parameter = parameter;
        this.lowerBound = lowerBound;
        this.upperBound = upperBound;
        this.offset = 0;
    }

    public abstract double[] getGradientLogDensity(double[] gradient);

    /**
     * Evaluate this gradient using an explicit source offset for this call.
     * This avoids relying on shared mutable offset state when the same
     * provider instance is reused across wrappers (e.g. path/source-destination
     * compositions).
     */
    public double[] getGradientLogDensity(double[] gradient, int sourceOffset) {
        final int previousOffset = this.offset;
        this.offset = sourceOffset;
        try {
            return getGradientLogDensity(gradient);
        } finally {
            this.offset = previousOffset;
        }
    }

    public abstract Parameter getRawParameter();

    public void setOffset(int offset) {
        this.offset = offset;
    }

    protected int getSourceOffset() {
        return offset;
    }

    protected double[] extractSourceGradient(double[] fullGradient, int requiredDimension) {
        final double[] result = new double[requiredDimension];
        if (fullGradient == null) {
            return result;
        }
        if (offset < 0) {
            throw new IllegalStateException("Negative source gradient offset: " + offset);
        }
        if (offset + requiredDimension > fullGradient.length) {
            throw new IllegalStateException(
                    "Source gradient slice out of bounds: offset=" + offset +
                            ", requiredDimension=" + requiredDimension +
                            ", fullGradientLength=" + fullGradient.length + "."
            );
        }
        System.arraycopy(fullGradient, offset, result, 0, requiredDimension);
        return result;
    }

    public abstract ContinuousTraitGradientForBranch.ContinuousProcessParameterGradient.DerivationParameter getDerivationParameter();

    public ContinuousTraitGradientForBranch.ProcessGradientTarget getProcessGradientTarget() {
        return ContinuousTraitGradientForBranch.targetForDerivationParameter(getDerivationParameter());
    }

    @Override
    public Likelihood getLikelihood() {
        return likelihood;
    }

//    @Override
//    public Parameter getParameter() {
//        return parameter;
//    }

    protected Parameter getNumericalParameter() {
        return getParameter();
    }

    @Override
    public String getReport() {
        return GradientWrtParameterProvider.getReportAndCheckForError(this,
                lowerBound, upperBound, TOLERANCE, REPORT_SMALL_COMPONENT_THRESHOLD)
                + "\n" + checkNumeric(getGradientLogDensity());
    }

//    @Override
//    public String getReport() {
//        return checkNumeric(getGradientLogDensity());
//    }

    String getReportString(double[] analytic, double[] numeric) {

        return getClass().getCanonicalName() + "\n" +
                "analytic: " + new Vector(analytic) +
                "\n" +
                "numeric: " + new Vector(numeric) +
                "\n";
    }

    String getReportString(double[] analytic, double[] numeric, double[] numericTrans) {

        return getClass().getCanonicalName() + "\n" +
                "analytic: " + new Vector(analytic) +
                "\n" +
                "numeric (no Cholesky): " + new Vector(numeric) +
                "\n" +
                "numeric (with Cholesky): " + new Vector(numericTrans) +
                "\n";
    }

    MultivariateFunction getNumeric() {

        return new MultivariateFunction() {

            @Override
            public double evaluate(double[] argument) {

                for (int i = 0; i < argument.length; ++i) {
                    getNumericalParameter().setParameterValue(i, argument[i]);
                }

                likelihood.makeDirty();
                System.err.println("likelihood in numeric:" + likelihood.getLogLikelihood());
                return likelihood.getLogLikelihood();
            }

            @Override
            public int getNumArguments() {
                return getDimension();
            }

            @Override
            public double getLowerBound(int n) {
                return lowerBound;
            }

            @Override
            public double getUpperBound(int n) {
                return upperBound;
            }
        };
    }

    String checkNumeric(double[] analytic) {

        System.err.println("Numeric at: \n" + new Vector(getNumericalParameter().getParameterValues()));

        double[] storedValues = getNumericalParameter().getParameterValues();
        double[] testGradient = NumericalDerivative.gradient(getNumeric(), storedValues);
        for (int i = 0; i < storedValues.length; ++i) {
            getNumericalParameter().setParameterValue(i, storedValues[i]);
        }

        return getReportString(analytic, testGradient);
    }

    double[] numericGradientRespectingBounds(double lower, double upper) {
        final Parameter parameter = getNumericalParameter();
        final double[] storedValues = parameter.getParameterValues();
        final double[] testGradient = new double[storedValues.length];

        for (int i = 0; i < storedValues.length; ++i) {
            final double x = storedValues[i];
            final double h = Math.pow(dr.math.MachineAccuracy.EPSILON, 0.333) * Math.max(Math.abs(x), 1.0);

            final double forward = x + h;
            final double backward = x - h;

            if (backward >= lower && forward <= upper) {
                storedValues[i] = forward;
                parameter.setParameterValue(i, storedValues[i]);
                likelihood.makeDirty();
                final double fxPlus = likelihood.getLogLikelihood();

                storedValues[i] = backward;
                parameter.setParameterValue(i, storedValues[i]);
                likelihood.makeDirty();
                final double fxMinus = likelihood.getLogLikelihood();

                testGradient[i] = (fxPlus - fxMinus) / (2.0 * h);
            } else if (forward <= upper) {
                final double localH = Math.min(h, upper - x);
                storedValues[i] = x;
                parameter.setParameterValue(i, storedValues[i]);
                likelihood.makeDirty();
                final double fx = likelihood.getLogLikelihood();

                storedValues[i] = x + localH;
                parameter.setParameterValue(i, storedValues[i]);
                likelihood.makeDirty();
                final double fxPlus = likelihood.getLogLikelihood();

                testGradient[i] = (fxPlus - fx) / localH;
            } else if (backward >= lower) {
                final double localH = Math.min(h, x - lower);
                storedValues[i] = x;
                parameter.setParameterValue(i, storedValues[i]);
                likelihood.makeDirty();
                final double fx = likelihood.getLogLikelihood();

                storedValues[i] = x - localH;
                parameter.setParameterValue(i, storedValues[i]);
                likelihood.makeDirty();
                final double fxMinus = likelihood.getLogLikelihood();

                testGradient[i] = (fx - fxMinus) / localH;
            } else {
                testGradient[i] = Double.NaN;
            }

            storedValues[i] = x;
            parameter.setParameterValue(i, x);
        }

        likelihood.makeDirty();
        return testGradient;
    }

    String checkNumericRespectingBounds(double[] analytic, double lower, double upper) {

        System.err.println("Numeric at: \n" + new Vector(getNumericalParameter().getParameterValues()));
        final double[] testGradient = numericGradientRespectingBounds(lower, upper);
        return "analytic: " + new Vector(analytic) +
                "\n" +
                "numeric : " + new Vector(testGradient) +
                "\n";
    }

    // Default branch specific class
    public static class ParameterDiffusionGradient extends AbstractDiffusionGradient implements Reportable {

        protected final int dim;
        protected final BranchSpecificGradient branchSpecificGradient;

        private final Parameter parameter;
        private final Parameter rawParameter;

        ParameterDiffusionGradient(BranchSpecificGradient branchSpecificGradient,
                                   Likelihood likelihood,
                                   Parameter parameter,
                                   Parameter rawParameter,
                                   double upperBound,
                                   double lowerBound) {
            super(likelihood, upperBound, lowerBound);

            this.parameter = parameter;
            this.rawParameter = rawParameter;

            this.branchSpecificGradient = branchSpecificGradient;
            this.dim = parameter.getDimension();

        }

        @Override
        public Parameter getParameter() {
            return parameter;
        }

        @Override
        public int getDimension() {
            return dim;
        }

        @Override
        public Parameter getRawParameter() {
            return rawParameter;
        }

        @Override
        public double[] getGradientLogDensity() {
            if (Boolean.getBoolean("beast.debug.parameterDiffusion.refreshLikelihoodBeforeAnalytic")) {
                getLikelihood().makeDirty();
                getLikelihood().getLogLikelihood();
            }
            double[] gradient = getFullBranchGradient();
            final double[] extracted = getGradientLogDensity(gradient);
            maybeEmitParameterDiffusionSliceDebug(gradient, extracted);
            maybeEmitRecomputeConsistencyDebug(gradient, extracted);
            return extracted;
        }

        protected double[] getFullBranchGradient() {
            if (Boolean.getBoolean("beast.gradientHook.parameterDiffusion.useFreshSnapshot")) {
                return branchSpecificGradient.getGradientLogDensityFromFreshLikelihoodSnapshot();
            }
            return branchSpecificGradient.getGradientLogDensity();
        }

        public double[] getGradientLogDensity(double[] gradient) {
            return extractSourceGradient(gradient, dim);
        }

        private void maybeEmitParameterDiffusionSliceDebug(final double[] fullAnalyticGradient,
                                                           final double[] extractedAnalyticSlice) {
            if (!Boolean.getBoolean("beast.debug.parameterDiffusionSlice")) {
                return;
            }
            final String rawName = rawParameter == null ? "<null>" : rawParameter.getParameterName();
            final String filter = System.getProperty(
                    "beast.debug.parameterDiffusionSlice.onlyParameterContains",
                    "OrthogonalBlockDiagonalPolarStableMatrixParameter.native"
            );
            if (filter != null && !filter.isEmpty() && !rawName.contains(filter)) {
                return;
            }
            final long call = DEBUG_PARAMETER_DIFFUSION_CALL_COUNTER.incrementAndGet();
            final int stride = Integer.getInteger("beast.debug.parameterDiffusionSliceStride", 1);
            if (stride <= 0 || (call % stride != 0)) {
                return;
            }

            final int[] focus = parseFocusIndices();
            final double[] fullNumericGradient = branchSpecificGradient.getNumericalGradientLogDensityDebug();
            final double[] extractedNumericSlice = extractSourceGradient(fullNumericGradient, dim);

            final StringBuilder sb = new StringBuilder();
            sb.append("parameterDiffusionSliceDebug call=").append(call)
                    .append(" rawParameter=").append(rawName)
                    .append(" sourceOffset=").append(getSourceOffset())
                    .append(" dim=").append(dim)
                    .append('\n');
            sb.append("\tfullAnalytic");
            appendIndexValues(sb, fullAnalyticGradient, focus);
            sb.append('\n');
            sb.append("\tfullNumeric ");
            appendIndexValues(sb, fullNumericGradient, focus);
            sb.append('\n');
            sb.append("\tfullAbsDiff");
            appendIndexAbsDiff(sb, fullAnalyticGradient, fullNumericGradient, focus);
            sb.append('\n');
            sb.append("\textractedAnalytic=").append(new Vector(extractedAnalyticSlice)).append('\n');
            sb.append("\textractedNumeric =").append(new Vector(extractedNumericSlice)).append('\n');
            sb.append("\textractedAbsDiff=").append(new Vector(absDiff(extractedAnalyticSlice, extractedNumericSlice))).append('\n');
            System.err.print(sb.toString());
        }

        /**
         * Single-state replay diagnostic used by HMC mismatch debugging.
         *
         * Compares, at the current parameter state:
         * 1) full branch analytic vs full branch numeric
         * 2) extracted native-block analytic vs extracted native-block numeric
         * and reports repeatability across multiple analytic evaluations.
         */
        public String getSingleStateReplayDebug() {
            final int repeats = Integer.getInteger("beast.debug.parameterDiffusionReplay.repeats", 3);
            final int[] focus = parseFocusIndices();

            final double[] fullAnalytic0 = branchSpecificGradient.getGradientLogDensity();
            final double[] fullNumeric = branchSpecificGradient.getNumericalGradientLogDensityDebug();
            final double[] extractedAnalytic0 = getGradientLogDensity(fullAnalytic0);
            final double[] extractedNumeric = getGradientLogDensity(fullNumeric);
            final String dumpPath = System.getProperty("beast.debug.parameterDiffusionReplay.dumpPath");
            final String dumpOnlyParameterContains = System.getProperty(
                    "beast.debug.parameterDiffusionReplay.dumpOnlyParameterContains", "");
            final String rawName = rawParameter == null ? "<null>" : rawParameter.getParameterName();
            boolean dumpAttempted = false;
            boolean dumpWritten = false;
            String dumpError = null;
            final boolean allowDumpForThisParameter =
                    dumpOnlyParameterContains.isEmpty() || rawName.contains(dumpOnlyParameterContains);
            if (allowDumpForThisParameter && dumpPath != null && !dumpPath.isEmpty()) {
                dumpAttempted = true;
                try {
                    branchSpecificGradient.dumpCurrentBranchStateForPythonDebug(dumpPath, fullAnalytic0, fullNumeric);
                    dumpWritten = true;
                } catch (RuntimeException ex) {
                    dumpError = ex.getMessage();
                }
            }

            double maxRepeatAbsDiff = 0.0;
            int maxRepeatIdx = -1;
            for (int r = 1; r < repeats; ++r) {
                final double[] fullAnalyticR = branchSpecificGradient.getGradientLogDensity();
                final double[] extractedAnalyticR = getGradientLogDensity(fullAnalyticR);
                for (int i = 0; i < extractedAnalytic0.length && i < extractedAnalyticR.length; ++i) {
                    final double d = Math.abs(extractedAnalytic0[i] - extractedAnalyticR[i]);
                    if (d > maxRepeatAbsDiff) {
                        maxRepeatAbsDiff = d;
                        maxRepeatIdx = i;
                    }
                }
            }

            final StringBuilder sb = new StringBuilder();
            sb.append("parameterDiffusionSingleStateReplay")
                    .append(" rawParameter=")
                    .append(rawParameter == null ? "<null>" : rawParameter.getParameterName())
                    .append(" sourceOffset=").append(getSourceOffset())
                    .append(" dim=").append(dim)
                    .append('\n');
            sb.append("\tfullAnalytic");
            appendIndexValues(sb, fullAnalytic0, focus);
            sb.append('\n');
            sb.append("\tfullNumeric ");
            appendIndexValues(sb, fullNumeric, focus);
            sb.append('\n');
            sb.append("\tfullAbsDiff");
            appendIndexAbsDiff(sb, fullAnalytic0, fullNumeric, focus);
            sb.append('\n');
            sb.append("\textractedAnalytic=").append(new Vector(extractedAnalytic0)).append('\n');
            sb.append("\textractedNumeric =").append(new Vector(extractedNumeric)).append('\n');
            sb.append("\textractedAbsDiff=").append(new Vector(absDiff(extractedAnalytic0, extractedNumeric))).append('\n');
            sb.append("\trepeatability maxAbsDiff=").append(maxRepeatAbsDiff)
                    .append(" maxIdx=").append(maxRepeatIdx)
                    .append(" repeats=").append(repeats)
                    .append('\n');
            if (dumpAttempted) {
                sb.append("\tdumpPath=").append(dumpPath)
                        .append(" dumpWritten=").append(dumpWritten);
                if (dumpError != null) {
                    sb.append(" dumpError=").append(dumpError.replace('\n', ' '));
                }
                sb.append('\n');
            }
            if (Boolean.getBoolean("beast.debug.parameterDiffusionReplay.comparePassStats")) {
                sb.append('\t')
                        .append(branchSpecificGradient.getConsecutivePassStatisticsDiffReport().replace("\n", "\n\t"))
                        .append('\n');
            }
            return sb.toString();
        }

        private int[] parseFocusIndices() {
            final String raw = System.getProperty("beast.debug.parameterDiffusionSlice.indices", "6,7,8");
            final String[] tokens = raw.split(",");
            final int[] out = new int[tokens.length];
            for (int i = 0; i < tokens.length; ++i) {
                try {
                    out[i] = Integer.parseInt(tokens[i].trim());
                } catch (NumberFormatException ignored) {
                    out[i] = -1;
                }
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

        private void appendIndexAbsDiff(final StringBuilder sb,
                                        final double[] left,
                                        final double[] right,
                                        final int[] indices) {
            for (int idx : indices) {
                if (idx >= 0 && idx < left.length && idx < right.length) {
                    sb.append(" i").append(idx).append('=').append(Math.abs(left[idx] - right[idx]));
                } else {
                    sb.append(" i").append(idx).append("=<oob>");
                }
            }
        }

        private double[] absDiff(final double[] a, final double[] b) {
            final int length = Math.min(a.length, b.length);
            final double[] out = new double[length];
            for (int i = 0; i < length; ++i) {
                out[i] = Math.abs(a[i] - b[i]);
            }
            return out;
        }

        private void maybeEmitRecomputeConsistencyDebug(final double[] baselineFullGradient,
                                                        final double[] baselineExtractedSlice) {
            if (!Boolean.getBoolean("beast.debug.parameterDiffusionRecomputeConsistency")) {
                return;
            }
            final String rawName = rawParameter == null ? "<null>" : rawParameter.getParameterName();
            final String filter = System.getProperty(
                    "beast.debug.parameterDiffusionRecomputeConsistency.onlyParameterContains",
                    ""
            );
            if (filter != null && !filter.isEmpty() && !rawName.contains(filter)) {
                return;
            }

            final long call = DEBUG_PARAMETER_DIFFUSION_CALL_COUNTER.incrementAndGet();
            final int stride = Integer.getInteger("beast.debug.parameterDiffusionRecomputeConsistencyStride", 1);
            if (stride <= 0 || (call % stride != 0)) {
                return;
            }

            getLikelihood().makeDirty();
            if (Boolean.getBoolean("beast.debug.parameterDiffusionRecomputeConsistencyForceLikelihood")) {
                getLikelihood().getLogLikelihood();
            }
            final double[] recomputedFullGradient = getFullBranchGradient();
            final double[] recomputedExtractedSlice = getGradientLogDensity(recomputedFullGradient);

            double maxAbsFull = 0.0;
            int maxIdxFull = -1;
            final int nFull = Math.min(baselineFullGradient.length, recomputedFullGradient.length);
            for (int i = 0; i < nFull; ++i) {
                final double d = Math.abs(baselineFullGradient[i] - recomputedFullGradient[i]);
                if (d > maxAbsFull) {
                    maxAbsFull = d;
                    maxIdxFull = i;
                }
            }

            double maxAbsSlice = 0.0;
            int maxIdxSlice = -1;
            final int nSlice = Math.min(baselineExtractedSlice.length, recomputedExtractedSlice.length);
            for (int i = 0; i < nSlice; ++i) {
                final double d = Math.abs(baselineExtractedSlice[i] - recomputedExtractedSlice[i]);
                if (d > maxAbsSlice) {
                    maxAbsSlice = d;
                    maxIdxSlice = i;
                }
            }

            System.err.println("parameterDiffusionRecomputeConsistency call=" + call
                    + " rawParameter=" + rawName
                    + " sourceOffset=" + getSourceOffset()
                    + " maxAbsFull=" + maxAbsFull
                    + " maxIdxFull=" + maxIdxFull
                    + " maxAbsSlice=" + maxAbsSlice
                    + " maxIdxSlice=" + maxIdxSlice);
        }

        @Override
        public ContinuousTraitGradientForBranch.ContinuousProcessParameterGradient.DerivationParameter getDerivationParameter() {
            List<ContinuousTraitGradientForBranch.ContinuousProcessParameterGradient.DerivationParameter> derList = branchSpecificGradient.getDerivationParameter();
            assert derList.size() == 1;
            return derList.get(0);
        }

        @Override
        public String getReport() {
            return "Gradient." + rawParameter.getParameterName() + "\n" +
                    super.getReport();
        }

        public static ParameterDiffusionGradient createDriftGradient(BranchSpecificGradient branchSpecificGradient,
                                                                     Likelihood likelihood,
                                                                     Parameter drift) {
            return new ParameterDiffusionGradient(
                    branchSpecificGradient, likelihood,
                    drift, drift,
                    Double.POSITIVE_INFINITY, Double.NEGATIVE_INFINITY);
        }

        public static ParameterDiffusionGradient createDiagonalAttenuationGradient(BranchSpecificGradient branchSpecificGradient,
                                                                                   Likelihood likelihood,
                                                                                   MatrixParameterInterface attenuation) {
            assert (attenuation instanceof DiagonalMatrix)
                    : "DiagonalAttenuationGradient can only be applied to a DiagonalMatrix.";

            return new ParameterDiffusionGradient(
                    branchSpecificGradient, likelihood,
                    ((DiagonalMatrix) attenuation).getDiagonalParameter(),
                    (DiagonalMatrix) attenuation,
                    Double.POSITIVE_INFINITY, 1.0e-12) {
                @Override
                public String getReport() {
                    final double[] analytic = getGradientLogDensity();
                    System.err.println("Numeric at: \n" + new Vector(getNumericalParameter().getParameterValues()));
                    final double[] numeric = numericGradientRespectingBounds(1.0e-12, Double.POSITIVE_INFINITY);
                    return "Gradient." + getRawParameter().getParameterName() + "\n" +
                            "Gradient\n" +
                            "analytic: " + new Vector(analytic) + "\n" +
                            "numeric : " + new Vector(numeric) + "\n\n" +
                            AbstractDiffusionGradient.ParameterDiffusionGradient.class.getCanonicalName() + "\n" +
                            "analytic: " + new Vector(analytic) + "\n" +
                            "numeric: " + new Vector(numeric) + "\n";
                }
            };
        }

        public static ParameterDiffusionGradient createSelectionParameterGradient(final BranchSpecificGradient branchSpecificGradient,
                                                                                 final Likelihood likelihood,
                                                                                 final Parameter parameter,
                                                                                 final AbstractBlockDiagonalTwoByTwoMatrixParameter rawParameter) {
            final boolean useFreshSnapshotHook =
                    Boolean.getBoolean("beast.gradientHook.selection.useFreshSnapshot");

            if (useFreshSnapshotHook) {
                return new ParameterDiffusionGradient(
                        branchSpecificGradient,
                        likelihood,
                        parameter,
                        rawParameter.getParameter(),
                        Double.POSITIVE_INFINITY,
                        Double.NEGATIVE_INFINITY) {
                    @Override
                    protected double[] getFullBranchGradient() {
                        return branchSpecificGradient.getGradientLogDensityFromFreshLikelihoodSnapshot();
                    }

                    @Override
                    public String getReport() {
                        return "Gradient." + getRawParameter().getParameterName() + " (freshSnapshotHook)\n" +
                                super.getReport();
                    }
                };
            }

            return new ParameterDiffusionGradient(
                    branchSpecificGradient,
                    likelihood,
                    parameter,
                    rawParameter.getParameter(),
                    Double.POSITIVE_INFINITY,
                    Double.NEGATIVE_INFINITY);
        }
    }
}
