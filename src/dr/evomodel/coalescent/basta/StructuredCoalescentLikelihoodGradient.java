/*
 * StructuredCoalescentLikelihoodGradient.java
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

import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.branchratemodel.DifferentiableBranchRates;
import dr.evomodel.coalescent.AbstractStructuredCoalescentLikelihood;
import dr.evomodel.coalescent.mascot.MascotLikelihood;
import dr.evomodel.substmodel.ComplexSubstitutionModelGradientSupport;
import dr.inference.hmc.GradientWrtParameterProvider;
import dr.evomodel.substmodel.SubstitutionModel;
import dr.inference.loggers.LogColumn;
import dr.inference.loggers.Loggable;
import dr.inference.model.*;
import dr.xml.Reportable;

/**
 * HMC-facing gradient provider shared by both structured-coalescent engines:
 * {@link BastaLikelihood} (matrix-exponential, {@code wrtParameter="migrationRate"}
 * or {@code "populationSize"}) and {@link MascotLikelihood} (RK4/adjoint ODE,
 * additionally {@code "clockRate"}). The two engines compute their gradients
 * in genuinely different intermediate representations -- BASTA's {@link
 * BastaLikelihoodDelegate#calculateGradient} returns a raw, un-chain-ruled
 * derivative over the {@code stateCount x stateCount} instantaneous-rate
 * matrix (or the {@code stateCount}-dimensioned population-size block), which
 * {@link WrtParameter#chainRule} then reduces to the free parameter; MASCOT's
 * {@code MascotCore} adjoint pass instead returns a gradient already expressed
 * directly in theta-space (raw migration rates, log population sizes), fully
 * reduced by {@link MascotLikelihood} itself ({@code
 * MascotDynamics#writeMigrationGradient}/{@code writePopSizeGradient}) before
 * this class ever sees it -- so {@code AbstractStructuredCoalescentLikelihood#
 * getGradientLogDensity(StructuredCoalescentLikelihoodGradient)} is the shared
 * template-method seam: {@link BastaLikelihood}'s override does the
 * intermediate-gradient-plus-{@link WrtParameter#chainRule} dance,
 * {@link MascotLikelihood}'s override just dispatches {@link #getType()} to
 * its own already-reduced per-part getters. This class itself only owns what
 * is common to both: the free {@link Parameter} HMC operates on, and (for
 * BASTA only) the {@link #chainRule}/{@link #requiresTransitionMatrices}/
 * {@link #getIntermediateGradientDimension} hooks {@link BastaLikelihoodDelegate}
 * needs.
 */
public class StructuredCoalescentLikelihoodGradient implements
        GradientWrtParameterProvider, ModelListener, Reportable, Loggable {

    private final AbstractStructuredCoalescentLikelihood structuredCoalescentLikelihood;
    private final SubstitutionModel substitutionModel;
    private final WrtParameter wrtParameter;

    private final Parameter parameter;

    private final int stateCount;

    private double numericGradientStepSize = GradientWrtParameterProvider.super.getNumericGradientStepSize();

    public StructuredCoalescentLikelihoodGradient(AbstractStructuredCoalescentLikelihood structuredCoalescentLikelihood,
                                                  SubstitutionModel substitutionModel,
                                                  WrtParameter wrtParameter) {
        this.structuredCoalescentLikelihood = structuredCoalescentLikelihood;
        this.substitutionModel = substitutionModel;
        this.wrtParameter = wrtParameter;

        this.parameter = wrtParameter.getParameter(structuredCoalescentLikelihood, substitutionModel);

        this.stateCount = structuredCoalescentLikelihood.getStateCount();
    }

    @Override
    public Likelihood getLikelihood() {
        return structuredCoalescentLikelihood;
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
        return structuredCoalescentLikelihood.getGradientLogDensity(this);
    }

    @Override
    public double getNumericGradientStepSize() {
        return numericGradientStepSize;
    }

    @Override
    public void setNumericGradientStepSize(double ratio) {
        this.numericGradientStepSize = ratio;
    }

    /** BASTA-only: reduces {@link BastaLikelihoodDelegate#calculateGradient}'s raw intermediate gradient to {@link #parameter}. */
    double[] chainRule(double[] gradient) {
        return wrtParameter.chainRule(gradient, structuredCoalescentLikelihood, substitutionModel);
    }

    /** BASTA-only. */
    boolean requiresTransitionMatrices() {
        return wrtParameter.requiresTransitionMatrices();
    }

    public WrtParameter getType() { return wrtParameter; }

    /** BASTA-only. */
    public int getIntermediateGradientDimension() {
        return wrtParameter.getIntermediateGradientDimension(stateCount);
    }

    @Override
    public LogColumn[] getColumns() {
        return new LogColumn[0];
    }

    @Override
    public void modelChangedEvent(Model model, Object object, int index) {

    }

    @Override
    public void modelRestored(Model model) {

    }


    @Override
    public String getReport() {

        StringBuilder sb = new StringBuilder();

        String message = GradientWrtParameterProvider.getReportAndCheckForError(this, 0.0, Double.POSITIVE_INFINITY, 10.0);
        sb.append(message);


        return  sb.toString();
    }

    public enum WrtParameter {
        MIGRATION_RATE("migrationRate") {
            @Override
            Parameter getParameter(AbstractStructuredCoalescentLikelihood structuredCoalescentLikelihood, SubstitutionModel substitutionModel) {
                if (structuredCoalescentLikelihood instanceof MascotLikelihood) {
                    MascotLikelihood mascotLikelihood = (MascotLikelihood) structuredCoalescentLikelihood;
                    String error = mascotLikelihood.getMigrationGradientCompatibilityError();
                    if (error != null) {
                        throw new IllegalArgumentException(error);
                    }
                    return mascotLikelihood.getMigrationRates();
                }
                String error = ComplexSubstitutionModelGradientSupport.getCompatibilityError(
                        substitutionModel, "structuredCoalescentLikelihoodGradient wrtParameter=\"migrationRate\"");
                if (error != null) {
                    throw new IllegalArgumentException(error);
                }
                return ComplexSubstitutionModelGradientSupport.getRatesParameter(substitutionModel);
            }

            @Override
            double[] chainRule(double[] gradient, AbstractStructuredCoalescentLikelihood structuredCoalescentLikelihood,
                               SubstitutionModel substitutionModel) {
                // BASTA-only: MascotLikelihood#getGradientLogDensity(StructuredCoalescentLikelihoodGradient)
                // never calls this -- MascotDynamics#writeMigrationGradient already returns a
                // fully-reduced gradient wrt the free migration-rate parameter.
                Parameter frequencies = substitutionModel.getFrequencyModel().getFrequencyParameter();
                final int dim = frequencies.getDimension();

                double[] chainedGradient = new double[dim * (dim - 1)];

                int k = 0;
                for (int i = 0; i < dim; ++i) {
                    for (int j = i + 1; j < dim; ++j) {
                        chainedGradient[k] = (gradient[i * dim + j] - gradient[i * dim + i]) * frequencies.getParameterValue(j);
                        ++k;
                    }
                }

                for (int j = 0; j < dim; ++j) {
                    for (int i = j + 1; i < dim; ++i) {
                        chainedGradient[k] = (gradient[i * dim + j] - gradient[i * dim + i]) * frequencies.getParameterValue(j);
                        ++k;
                    }
                }
                return chainedGradient;
            }

            @Override
            int getIntermediateGradientDimension(int stateCount) {
                return stateCount * stateCount;
            }

            @Override
            boolean requiresTransitionMatrices() {
                return true;
            }
        },

        POPULATION_SIZE("populationSize") {
            @Override
            Parameter getParameter(AbstractStructuredCoalescentLikelihood structuredCoalescentLikelihood, SubstitutionModel substitutionModel) {
                if (structuredCoalescentLikelihood instanceof MascotLikelihood) {
                    return ((MascotLikelihood) structuredCoalescentLikelihood).getPopSizes();
                }
                return bastaPopulationSizeParameter(structuredCoalescentLikelihood);
            }

            @Override
            double[] chainRule(double[] gradient, AbstractStructuredCoalescentLikelihood structuredCoalescentLikelihood,
                               SubstitutionModel substitutionModel) {
                // BASTA-only: see the MIGRATION_RATE note above; MASCOT's
                // MascotDynamics#writePopSizeGradient already applies the
                // natural-scale chain rule before this class sees it.
                Parameter popSizes = bastaPopulationSizeParameter(structuredCoalescentLikelihood);
                final int dim = popSizes.getDimension();

                for (int i = 0; i < dim; ++i) {
                    double popSize = popSizes.getParameterValue(i);
                    gradient[i] /= -(popSize * popSize);
                }

                return gradient;
            }

            @Override
            int getIntermediateGradientDimension(int stateCount) {
                return stateCount;
            }

            @Override
            boolean requiresTransitionMatrices() {
                return false;
            }
        },

        CLOCK_RATE("clockRate") {
            @Override
            Parameter getParameter(AbstractStructuredCoalescentLikelihood structuredCoalescentLikelihood, SubstitutionModel substitutionModel) {
                if (!(structuredCoalescentLikelihood instanceof MascotLikelihood)) {
                    throw new IllegalArgumentException(
                            "structuredCoalescentLikelihoodGradient wrtParameter=\"clockRate\" is only supported for a MASCOT likelihood");
                }
                BranchRateModel branchRateModel = structuredCoalescentLikelihood.getBranchRateModel();
                if (!(branchRateModel instanceof DifferentiableBranchRates)) {
                    throw new IllegalArgumentException(
                            "structuredCoalescentLikelihoodGradient wrtParameter=\"clockRate\" requires a branchRateModel " +
                                    "that implements DifferentiableBranchRates, got: " +
                                    (branchRateModel == null ? "null (no branchRateModel supplied)" :
                                            branchRateModel.getClass().getName()));
                }
                return ((DifferentiableBranchRates) branchRateModel).getRateParameter();
            }

            @Override
            double[] chainRule(double[] gradient, AbstractStructuredCoalescentLikelihood structuredCoalescentLikelihood,
                               SubstitutionModel substitutionModel) {
                // Unreachable: getParameter() above throws for any non-MASCOT likelihood, and
                // MascotLikelihood#getGradientLogDensity(StructuredCoalescentLikelihoodGradient) returns
                // the clock-rate gradient directly (see MascotLikelihood#getClockRateGradientLogDensity())
                // rather than going through this BASTA-only intermediate-gradient-plus-chainRule path.
                throw new UnsupportedOperationException("wrtParameter=\"clockRate\" does not use chainRule()");
            }

            @Override
            int getIntermediateGradientDimension(int stateCount) {
                throw new UnsupportedOperationException("wrtParameter=\"clockRate\" has no intermediate gradient");
            }

            @Override
            boolean requiresTransitionMatrices() {
                throw new UnsupportedOperationException("wrtParameter=\"clockRate\" has no transition-matrix requirement");
            }
        };

        WrtParameter(String name) {
            this.name = name;
        }

        abstract Parameter getParameter(AbstractStructuredCoalescentLikelihood structuredCoalescentLikelihood, SubstitutionModel substitutionModel);

        abstract double[] chainRule(double[] gradient, AbstractStructuredCoalescentLikelihood structuredCoalescentLikelihood,
                                    SubstitutionModel substitutionModel);

        abstract int getIntermediateGradientDimension(int stateCount);

        abstract boolean requiresTransitionMatrices();

        private final String name;

        public static WrtParameter factory(String match) {
            for (WrtParameter type : WrtParameter.values()) {
                if (match.equalsIgnoreCase(type.name)) {
                    return type;
                }
            }
            return null;
        }

        /** Shared by {@code POPULATION_SIZE}'s {@code getParameter}/{@code chainRule} for a BASTA likelihood. */
        private static Parameter bastaPopulationSizeParameter(AbstractStructuredCoalescentLikelihood structuredCoalescentLikelihood) {
            AbstractPopulationSizeModel popModel = structuredCoalescentLikelihood.getPopulationSizeModel();
            if (popModel instanceof ConstantPopulationSizeModel) {
                return ((ConstantPopulationSizeModel) popModel).getPopulationSizeParameter();
            } else if (popModel instanceof ExponentialGrowthPopulationSizeModel) {
                return ((ExponentialGrowthPopulationSizeModel) popModel).getPopulationSizeParameter();
            } else if (popModel instanceof AnchoredExponentialGrowthPopulationSizeModel) {
                return ((AnchoredExponentialGrowthPopulationSizeModel) popModel).getLogPopSizesParameter();
            } else if (popModel instanceof SharedAncestralExponentialGrowthPopulationSizeModel) {
                throw new UnsupportedOperationException("Analytic population-size gradient is not implemented for the " +
                        "shared-ancestral (divergence-anchored) exponential model; run without the population-size gradient.");
            }
            throw new RuntimeException("Unknown population size model type");
        }
    }
}
