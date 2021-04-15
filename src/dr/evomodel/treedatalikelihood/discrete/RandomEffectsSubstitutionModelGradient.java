/*
 * DiscreteTraitBranchRateGradient.java
 *
 * Copyright (c) 2002-2020 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.evomodel.treedatalikelihood.discrete;

import dr.evomodel.substmodel.OldGLMSubstitutionModel;
import dr.evomodel.treedatalikelihood.BeagleDataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.inference.distribution.GeneralizedLinearModel;
import dr.inference.model.CompoundParameter;
import dr.inference.model.DesignMatrix;
import dr.inference.model.Parameter;

import java.util.ArrayList;
import java.util.List;

/**
 * @author Marc A. Suchard
 */
public class RandomEffectsSubstitutionModelGradient extends GlmSubstitutionModelGradient {

    private final int[][] mapEffectToIndices;

    public RandomEffectsSubstitutionModelGradient(String traitName,
                                                  TreeDataLikelihood treeDataLikelihood,
                                                  BeagleDataLikelihoodDelegate likelihoodDelegate,
                                                  OldGLMSubstitutionModel substitutionModel) {
        super(traitName, treeDataLikelihood, likelihoodDelegate, substitutionModel);

        // Count random effects dimension
        int asymmetricCount = stateCount * (stateCount - 1);
        if (getDimension() == asymmetricCount) {
            // Asymmetric
            mapEffectToIndices = makeAsymmetricMap();
        } else if (getDimension() == asymmetricCount / 2) {
            // Symmetric
            throw new RuntimeException("Not yet implemented");
        } else {
            throw new IllegalArgumentException("Unable to determine random effects count");
        }
    }

//    Parameter makeCompoundParameter(GeneralizedLinearModel glm) {
//        CompoundParameter parameter = new CompoundParameter("random.effects");
//        for (int i = 0; i < glm.getNumberOfRandomEffects(); ++i) {
//            parameter.addParameter(glm.getRandomEffect(i));
//        }
//        return parameter;
//    }

    ParameterMap makeParameterMap(GeneralizedLinearModel glm) {

        return new ParameterMap() {
            public double[] getCovariateColumn(int i) {
                return null;
            }

            public Parameter getParameter() { return glm.getRandomEffect(0); }
        };
    }

    @Override
    double preProcessNormalization(double[] differentials, double[] generator,
                                   boolean normalize) {
        double total = 0.0;
        if (normalize) {
            for (int i = 0; i < stateCount; ++i) {
                for (int j = 0; j < stateCount; ++j) {
                    total += differentials[index(i, j)] * generator[index(i, j)];
                }
            }
        }
        return total;
    }

    @Override
    double processSingleGradientDimension(int k,
                                          double[] differentials, double[] generator, double[] pi,
                                          boolean normalize, double normalizationConstant) {

        double elementUpper = generator[indexIJ(k)];
        double total = (differentials[indexIJ(k)]  - differentials[indexII(k)]) * elementUpper;

        if (normalize) {
            total -= elementUpper * pi[indexI(k)] * normalizationConstant;
        }
        
        return total;
    }

    private int[][] makeAsymmetricMap() {
        int[][] map = new int[stateCount * (stateCount - 1)][];

        int k = 0;
        for (int i = 0; i < stateCount; ++i) {
            for (int j = i + 1; j < stateCount; ++j) {
                map[k++] = new int[]{i, j};
            }
        }

        for (int j = 0; j < stateCount; ++j) {
            for (int i = j + 1; i < stateCount; ++i) {
                map[k++] = new int[]{i, j};
            }
        }

        return map;
    }

    private int indexIJ(int k) {
        final int[] indices = mapEffectToIndices[k];
        return indices[0] * stateCount + indices[1];
    }

    private int indexII(int k) {
        final int[] indices = mapEffectToIndices[k];
        return indices[0] * stateCount + indices[0];
    }

    private int indexI(int k) {
        final int[] indices = mapEffectToIndices[k];
        return indices[0];
    }
}
