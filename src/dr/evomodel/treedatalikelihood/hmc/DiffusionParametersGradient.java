/*
 * DiffusionParametersGradient.java
 *
 * Copyright (c) 2002-2019 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.evomodel.treedatalikelihood.hmc;

import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.evomodel.treedatalikelihood.continuous.BranchSpecificGradient;
import dr.inference.hmc.CompoundGradient;
import dr.inference.hmc.GradientWrtParameterProvider;
import dr.inference.model.CompoundParameter;
import dr.inference.model.Likelihood;
import dr.inference.model.Parameter;
import dr.xml.Reportable;

public class DiffusionParametersGradient implements GradientWrtParameterProvider, Reportable {

    private final TreeDataLikelihood likelihood;
    private final int dim;
    private final BranchSpecificGradient branchSpecificGradient;
    private final CompoundParameter compoundParameter;
    private final CompoundGradient parametersGradients;

    public DiffusionParametersGradient(BranchSpecificGradient branchSpecificGradient, CompoundGradient parametersGradients) {

        this.branchSpecificGradient = branchSpecificGradient;
        this.likelihood = (TreeDataLikelihood) branchSpecificGradient.getLikelihood();

        compoundParameter = new CompoundParameter(null);
        dim = checkAndSetParametersGradients(parametersGradients, compoundParameter);
        this.parametersGradients = parametersGradients;

    }

    private int checkAndSetParametersGradients(CompoundGradient parametersGradients, CompoundParameter parameter) {
        int offset = 0;
        int dim = 0;
        int dimTrait = likelihood.getDataLikelihoodDelegate().getTraitDim();
        for (GradientWrtParameterProvider gradient : parametersGradients.getDerivativeList()) {
            assert gradient instanceof AbstractDiffusionGradient : "Gradients must all be instances of AbstractDiffusionGradient.";
            ((AbstractDiffusionGradient) gradient).setOffset(offset);
            parameter.addParameter(gradient.getParameter());
            offset += ((AbstractDiffusionGradient) gradient).getDerivationParameter().getDimension(dimTrait);
            dim += gradient.getDimension();
        }
        return dim;
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

        int offset = 0;
        for (GradientWrtParameterProvider gradientProvider : parametersGradients.getDerivativeList()) {
            System.arraycopy(((AbstractDiffusionGradient) gradientProvider).getGradientLogDensity(gradient), 0,
                    result, offset, gradientProvider.getDimension());
            offset += gradientProvider.getDimension();
        }
        return result;
    }

    @Override
    public String getReport() {
        return "diffusionGradient." + compoundParameter.getParameterName() + "\n" +
                GradientWrtParameterProvider.getReportAndCheckForError(this,
                        Double.NEGATIVE_INFINITY, Double.POSITIVE_INFINITY,
                        GradientWrtParameterProvider.TOLERANCE);
    }
}
