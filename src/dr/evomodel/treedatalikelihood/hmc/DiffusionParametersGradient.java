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
import java.util.Collections;
import java.util.List;

public class DiffusionParametersGradient implements GradientWrtParameterProvider, Reportable {

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
        double[] gradient = branchSpecificGradient.getGradientLogDensity();
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
        return result;
    }

    @Override
    public String getReport() {
        return "diffusionGradient." + compoundParameter.getParameterName() + "\n" +
                GradientWrtParameterProvider.getReportAndCheckForError(this,
                        Double.NEGATIVE_INFINITY, Double.POSITIVE_INFINITY,
                GradientWrtParameterProvider.TOLERANCE,
                GradientWrtParameterProvider.SMALL_NUMBER_THRESHOLD);
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

}
