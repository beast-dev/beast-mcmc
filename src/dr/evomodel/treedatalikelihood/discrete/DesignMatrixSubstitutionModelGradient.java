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
import dr.inference.model.DesignMatrix;
import dr.inference.model.MaskedParameter;
import dr.inference.model.Parameter;

/**
 * @author Marc A. Suchard
 */
public class DesignMatrixSubstitutionModelGradient extends GlmSubstitutionModelGradient {

    private final int whichEffect;
    private final MaskedParameter parameter;
    private final int[][] mapEffectToIndices;

    public DesignMatrixSubstitutionModelGradient(String traitName,
                                                 TreeDataLikelihood treeDataLikelihood,
                                                 BeagleDataLikelihoodDelegate likelihoodDelegate,
                                                 OldGLMSubstitutionModel substitutionModel,
                                                 DesignMatrix matrix,
                                                 MaskedParameter parameter) {
        super(traitName, treeDataLikelihood, likelihoodDelegate, substitutionModel);

        this.parameter = parameter;

        int dim = parameter.getUnmaskedParameter().getDimension();
        int asymmetricCount = stateCount * (stateCount - 1);
        if (dim == asymmetricCount) {
            // Asymmetric
            mapEffectToIndices = makeDesignMap();
        } else if (getDimension() == asymmetricCount / 2) {
            // Symmetric
            throw new RuntimeException("Not yet implemented");
        } else {
            throw new IllegalArgumentException("Unable to determine random design matrix count");
        }

        whichEffect = findDesignMatrix(matrix);

        if (glm.getFixedEffect(whichEffect).getDimension() != 1 ||
                (glm.getFixedEffectIndicator(whichEffect) != null &&
                glm.getFixedEffectIndicator(whichEffect).getDimension() != 1)) {
            throw new IllegalArgumentException("Not yet implemented");
        }
    }

    private int findDesignMatrix(DesignMatrix matrix) {
        for (int i = 0; i < glm.getNumberOfFixedEffects(); ++i) {
            if (glm.getDesignMatrix(i) == matrix) {
                return i;
            }
        }
        throw new IllegalArgumentException("Unable to find design matrix in GLM model");
    }

    Double getReportTolerance() {
        return 0.01; // TODO Change back to null
    }

    @Override
    public Parameter getParameter() {
        return parameter;
    }

    ParameterMap makeParameterMap(GeneralizedLinearModel glm) {
        return null;
    }

    @Override
    double preProcessNormalization(double[] differentials, double[] generator,
                                   boolean normalize) {
        // TODO Code duplication with RandomEffectsGlmSubstitutionModel
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


        double coefficient = glm.getFixedEffect(whichEffect).getParameterValue(0);
        Parameter indicator = glm.getFixedEffectIndicator(whichEffect);
        if (indicator != null) {
            coefficient *= indicator.getParameterValue(0);
        }

        double element = generator[indexIJ(k)] * coefficient;
        double total = (differentials[indexIJ(k)]  - differentials[indexII(k)]) * element;

        if (normalize) {
            total -= element * pi[indexI(k)] * normalizationConstant;
        }
        
        return total;
    }

    private int[][] makeDesignMap() {
        int[][] map = new int[parameter.getDimension()][];

        int k = 0, l = 0;
        for (int i = 0; i < stateCount; ++i) {
            for (int j = i + 1; j < stateCount; ++j) {
                if (parameter.getParameterMaskValue(l++) == 1.0) {
                    map[k++] = new int[]{i, j};
                }
            }
        }

        for (int j = 0; j < stateCount; ++j) {
            for (int i = j + 1; i < stateCount; ++i) {
                if (parameter.getParameterMaskValue(l++) == 1.0) {
                    map[k++] = new int[]{i, j};
                }
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
