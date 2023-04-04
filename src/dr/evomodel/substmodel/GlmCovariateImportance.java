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
import dr.inference.model.DesignMatrix;
import dr.inference.model.Likelihood;
import dr.inference.model.Parameter;
import dr.math.MathUtils;
import dr.xml.Reportable;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * @author Marc A. Suchard
 * @author Philippe Lemey
 */
public class GlmCovariateImportance implements Loggable, Reportable {

    private final Likelihood likelihood;
    private final GeneralizedLinearModel glm;

    private final int dim;
    private final List<Parameter> covariates;

    private double currentLogLikelihood = Double.NaN;
    private LogColumn[] columns;

    class IndexPair {
        int fixedEffectIndex;
        int parameterValueIndex;

        IndexPair(int effect, int index) {
            this.fixedEffectIndex = effect;
            this.parameterValueIndex = index;
        }
    }

    private final Map<Integer, IndexPair> parameterMap;

    public GlmCovariateImportance(Likelihood likelihood,
                                  OldGLMSubstitutionModel glmSubstitutionModel) {
        this.likelihood = likelihood;
        this.glm = glmSubstitutionModel.getGeneralizedLinearModel();
        covariates = new ArrayList<>();
        parameterMap = new HashMap<>();
        int dim = 0;
        for (int i = 0; i < glm.getNumberOfFixedEffects(); ++i) {
            DesignMatrix design = glm.getDesignMatrix(i);
            int d = design.getColumnDimension();
            for (int j = 0; j < d; ++j) {
                covariates.add(design.getParameter(j));
                parameterMap.put(dim + j, new IndexPair(i, j));
            }
            dim += d;
        }
        this.dim = dim;
    }

    private double getDeviance(int index) {
        if (index == 0) {
            if (!Double.isNaN(currentLogLikelihood)) {
                throw new IllegalStateException("GlmCovariateImportance computation is already in-process");
            }

            currentLogLikelihood = likelihood.getLogLikelihood();
        } else if (Double.isNaN(currentLogLikelihood)) {
            throw new IllegalStateException("GlmCovariateImportance computation is already in-process");
        }

        boolean doPermutation = true;
        if (index != dim) {

            IndexPair indexPair = parameterMap.get(index);
            double coefficient = glm.getFixedEffect(indexPair.fixedEffectIndex).getParameterValue(indexPair.parameterValueIndex);
            if (glm.getFixedEffectIndicator(indexPair.fixedEffectIndex) != null) {
                coefficient *= glm.getFixedEffectIndicator(indexPair.fixedEffectIndex).getParameterValue(indexPair.parameterValueIndex);
            }

            if (coefficient == 0.0) {
                doPermutation = false;
            }
        }

        final double deviance;
        if (doPermutation) {

            glm.storeModelState();
            permute(index);
            deviance = currentLogLikelihood - likelihood.getLogLikelihood();
            glm.restoreModelState();
            glm.fireModelChanged();
            
        } else {
            deviance = 0.0;
        }

        if (index == dim) {
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
        Parameter covariate = covariates.get(index);
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
                String name = "deviance" + (index < dim ? index + 1 : "All");
                columns[i + 1] = new NumberColumn(name) {
                    @Override
                    public double getDoubleValue() {
                        return getDeviance(index);
                    }
                };
            }
        }
        return columns;
    }

    @Override
    public String getReport() {

        StringBuilder sb = new StringBuilder();
        LogColumn[] columns = getColumns();
        for (int i = 0; i < columns.length; ++i) {
            String label = columns[i].getLabel();
            String value = columns[i].getFormatted();
            System.err.println(label + ": " + value);
            sb.append(label).append(": ").append(value).append("\n");
        }
        return sb.toString();
    }
}
