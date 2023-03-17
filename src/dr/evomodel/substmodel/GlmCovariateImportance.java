/*
 * GlmCovariateImportance.java
 *
 * Copyright (c) 2002-2023 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.evomodel.substmodel;

import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.inference.distribution.GeneralizedLinearModel;
import dr.inference.loggers.LogColumn;
import dr.inference.loggers.Loggable;
import dr.inference.loggers.NumberColumn;
import dr.inference.model.Likelihood;
import dr.inference.model.Parameter;
import dr.math.MathUtils;

/**
 * @author Marc A. Suchard
 * @author Philippe Lemey
 */
public class GlmCovariateImportance implements Loggable {

    private final Likelihood likelihood;
    private final GeneralizedLinearModel glm;
    private final int dim;

    private double currentLogLikelihood = Double.NaN;
    private LogColumn[] columns;

    public GlmCovariateImportance(TreeDataLikelihood likelihood,
                                  OldGLMSubstitutionModel glmSubstitutionModel) {
        this.likelihood = likelihood;
        this.glm = glmSubstitutionModel.getGeneralizedLinearModel();
        this.dim = glm.getNumberOfFixedEffects();
    }

    private double getDeviance(int index) {
        if (index == 0) {
            if (!Double.isNaN(currentLogLikelihood)) {
                throw new IllegalStateException("GlmCovariateImportance computation is already in-process");
            }

            currentLogLikelihood = likelihood.getLogLikelihood();
            glm.storeModelState();
        } else if (Double.isNaN(currentLogLikelihood)) {
            throw new IllegalStateException("GlmCovariateImportance computation is already in-process");
        }

        permute(index);
        double deviance = likelihood.getLogLikelihood() - currentLogLikelihood;

        if (index == dim) {
            glm.restoreModelState();
            likelihood.makeDirty(); // TODO may not be necessary
            double test = likelihood.getLogLikelihood();
            if (test != currentLogLikelihood) {
                throw new IllegalStateException("Unable to restore model after covariate permutation");
            }
            currentLogLikelihood = Double.NaN;
        }

        return deviance;
    }

    private void permute(int index) {
        if (index != dim) {
            permuteOneImpl(index);
        } else {
            for (int i = 0; i < dim; ++i) {
                permuteOneImpl(i);
            }
        }
    }

    private void permuteOneImpl(int index) {
        Parameter covariate = glm.getFixedEffect(index);
        double[] values = covariate.getParameterValues();
        shuffle(values);

        for (int i = 0; i < values.length; ++i) {
            covariate.setParameterValueQuietly(i, values[i]);
        }
        covariate.fireParameterChangedEvent();
    }

    private void shuffle(double[] array) {
        for (int i = array.length - 1; i > 0; i--) {
            int index = MathUtils.nextInt(i + 1);
            double temp = array[index];
            array[index] = array[i];
            array[i] = temp;
        }
    }

    @Override
    public LogColumn[] getColumns() {
        if (columns == null) {
            columns = new LogColumn[dim + 2];

            columns[0] = new NumberColumn("logLikelihood") {
                @Override
                public double getDoubleValue() {
                    return likelihood.getLogLikelihood();
                }
            };

            for (int i = 0; i < dim + 1; ++i) {
                final int index = i;
                columns[i + 1] = new NumberColumn("deviance" + index) {
                    @Override
                    public double getDoubleValue() {
                        return getDeviance(index);
                    }
                };
            }
        }
        return columns;
    }
}
