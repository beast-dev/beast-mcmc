/*
 * FullyConjugateTreeTipsPotentialDerivative.java
 *
 * Copyright (c) 2002-2017 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.evomodel.continuous.hmc;

import dr.evomodel.continuous.FullyConjugateMultivariateTraitLikelihood;
import dr.inference.hmc.GradientWrtParameterProvider;
import dr.inference.model.Likelihood;
import dr.inference.model.Parameter;

/**
 * @author Max Tolkoff
 */
public class FullyConjugateTreeTipsPotentialDerivative implements GradientWrtParameterProvider {

    private final FullyConjugateMultivariateTraitLikelihood treeLikelihood;
    private final Parameter traitParameter;

    public FullyConjugateTreeTipsPotentialDerivative(FullyConjugateMultivariateTraitLikelihood treeLikelihood){
        this.treeLikelihood = treeLikelihood;
        traitParameter = treeLikelihood.getTraitParameter();
    }

    @Override
    public Likelihood getLikelihood() {
        return treeLikelihood;
    }

    @Override
    public Parameter getParameter() {
        return traitParameter;
    }

    @Override
    public int getDimension() {
        return traitParameter.getDimension();
    }
    
    @Override
    public double[] getGradientLogDensity() {

        final int dimTraits = treeLikelihood.getDimTrait() * treeLikelihood.getNumData();
        final int ntaxa = traitParameter.getDimension() / dimTraits;

        final double[] derivative = new double[traitParameter.getDimension()];

        final double[][] allMeans = treeLikelihood.getConditionalMeans();
        final double[] allScalars = treeLikelihood.getPrecisionFactors();
        final double[][] precisionMatrix = treeLikelihood.getDiffusionModel().getPrecisionmatrix();

        for (int i = 0; i < ntaxa; ++i) {

            final double[] mean = allMeans[i];
            final double scale = allScalars[i];

            for (int j = 0; j < dimTraits; ++j) {

                double sum = 0.0;
                for (int k = 0; k < dimTraits; ++k) {
                    sum += (mean[k] - traitParameter.getParameterValue(i * dimTraits + k)) * scale * precisionMatrix[j][k];
                }
                derivative[i * dimTraits + j] = sum;
            }
        }

//        for (int i = 0; i < dimTraits; i++) { // This only works for IDENTITY matrices
//            for (int j = 0; j < ntaxa; j++) {
//                derivative[j * dimTraits + i] -= (traitParameter.getParameterValue(j * dimTraits + i) - mean[j][i]) * precfactor[j];
//                /* Sign change */
//            }
//
//        }

        return derivative;
    }
}
