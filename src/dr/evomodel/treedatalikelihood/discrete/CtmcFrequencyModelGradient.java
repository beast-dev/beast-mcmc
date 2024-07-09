/*
 * CtmcFrequencyModelGradient.java
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

package dr.evomodel.treedatalikelihood.discrete;

import dr.evomodel.substmodel.FrequencyModel;
import dr.evomodel.substmodel.GlmSubstitutionModel;
import dr.evomodel.substmodel.SubstitutionModel;
import dr.evomodel.treedatalikelihood.BeagleDataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.inference.loggers.LogColumn;
import dr.inference.model.Parameter;
import dr.util.Citation;

import java.util.List;

/**
 * @author Marc A. Suchard
 */

public class CtmcFrequencyModelGradient extends AbstractLogAdditiveSubstitutionModelGradient {

    private final FrequencyModel frequencyModel;

    public CtmcFrequencyModelGradient(String traitName,
                                      TreeDataLikelihood treeDataLikelihood,
                                      BeagleDataLikelihoodDelegate likelihoodDelegate,
                                      GlmSubstitutionModel substitutionModel) {
        super(traitName, treeDataLikelihood, likelihoodDelegate, substitutionModel,
                ApproximationMode.FIRST_ORDER);

        List<SubstitutionModel> substitutionModels = likelihoodDelegate.getBranchModel().getSubstitutionModels();
        this.frequencyModel = likelihoodDelegate.getBranchModel().getRootFrequencyModel();

        for (SubstitutionModel model : substitutionModels) {
            if (frequencyModel != model.getFrequencyModel()) {
                throw new RuntimeException("Not yet implemented");
            }
        }
    }

     @Override
    protected double preProcessNormalization(double[] differentials, double[] generator,
                                             boolean normalize) {
        double total = 0.0;
        if (normalize) {
            for (int i = 0; i < stateCount; ++i) {
                for (int j = 0; j < stateCount; ++j) {
                    final int ij = i * stateCount + j;
                    total += differentials[ij] * generator[ij];
                }
            }
        }
        return total;
    }

    @Override
    double processSingleGradientDimension(int j, double[] differentials, double[] generator, double[] pi,
                                          boolean normalize, double normalizationConstant) {
        // derivative wrt pi[j]
        double total = 0.0;

        for (int i = 0; i < stateCount; ++i) {
            final int ii = i * stateCount + i;
            final int ij = i * stateCount + j;
            total += (differentials[ij] - differentials[ii]) * generator[ij];
        }

        if (normalize) {
            for (int i = 0; i < stateCount; ++i) {
                final int ij = i * stateCount + j;
                total -= generator[ij] * pi[j] * normalizationConstant;
            }
        }

        return total;
    }

    @Override
    public Parameter getParameter() {
        return frequencyModel.getFrequencyParameter();
    }

    @Override
    public LogColumn[] getColumns() {
        throw new RuntimeException("Not yet implemented");
    }

    @Override
    public Citation.Category getCategory() {
        return Citation.Category.SUBSTITUTION_MODELS;
    }

    @Override
    public String getDescription() {
        return null; // TODO
    }

    @Override
    public List<Citation> getCitations() {
        // TODO Update
        return null;
    }
}
