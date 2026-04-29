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

public abstract class AbstractDiffusionGradient implements GradientWrtParameterProvider, Reportable {

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
                lowerBound, upperBound, TOLERANCE,
                GradientWrtParameterProvider.SMALL_NUMBER_THRESHOLD);
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
            double[] gradient = getFullBranchGradient();
            return getGradientLogDensity(gradient);
        }

        protected double[] getFullBranchGradient() {
            return branchSpecificGradient.getGradientLogDensity();
        }

        public double[] getGradientLogDensity(double[] gradient) {
            return extractSourceGradient(gradient, dim);
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
                    Double.POSITIVE_INFINITY, 1.0e-12);
        }

        public static ParameterDiffusionGradient createSelectionParameterGradient(final BranchSpecificGradient branchSpecificGradient,
                                                                                 final Likelihood likelihood,
                                                                                 final Parameter parameter,
                                                                                 final AbstractBlockDiagonalTwoByTwoMatrixParameter rawParameter) {
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
