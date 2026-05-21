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
import dr.evomodel.substmodel.BaseSubstitutionModel;
import dr.evomodel.substmodel.ComplexSubstitutionModel;
import dr.evomodel.substmodel.GeneralSubstitutionModel;
import dr.inference.hmc.GradientWrtParameterProvider;
import dr.evomodel.substmodel.SubstitutionModel;
import dr.inference.loggers.LogColumn;
import dr.inference.loggers.Loggable;
import dr.inference.model.*;
import dr.xml.Reportable;

public class StructuredCoalescentLikelihoodGradient implements
        GradientWrtParameterProvider, ModelListener, Reportable, Loggable {

    private final BastaLikelihood structuredCoalescentLikelihood;
    private final WrtParameter wrtParameter;
    private final SubstitutionModel substitutionModel;

    private final Parameter parameter;
    private final Parameter chainRuleDependent;

    private final int stateCount;

    public StructuredCoalescentLikelihoodGradient(BastaLikelihood BastaLikelihood,
                                                  SubstitutionModel substitutionModel,
                                                  WrtParameter wrtParameter) {
        this.structuredCoalescentLikelihood = BastaLikelihood;
        this.wrtParameter = wrtParameter;
        this.substitutionModel = substitutionModel;

        this.parameter = wrtParameter.getParameter(structuredCoalescentLikelihood, substitutionModel);
        this.chainRuleDependent = wrtParameter.getChainRuleDependent(structuredCoalescentLikelihood, substitutionModel);

        this.stateCount = structuredCoalescentLikelihood.getSubstitutionModel().getFrequencyModel().getFrequencyCount();
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

    double[] chainRule(double[] gradient) {
        return wrtParameter.chainRule(gradient, chainRuleDependent, substitutionModel);
    }

    boolean requiresTransitionMatrices() {
        return wrtParameter.requiresTransitionMatrices();
    }

    WrtParameter getType() { return wrtParameter; }

    public int getIntermediateGradientDimension() {
//        return structuredCoalescentLikelihood.
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
            Parameter getParameter(BastaLikelihood structuredCoalescentLikelihood, SubstitutionModel substitutionModel) {
                assert(substitutionModel instanceof GeneralSubstitutionModel);
                GeneralSubstitutionModel generalSubstitutionModel = (GeneralSubstitutionModel) substitutionModel;
                return generalSubstitutionModel.getRatesParameter();
            }

            @Override
            double[] chainRule(double[] gradient, Parameter parameter, SubstitutionModel substitutionModel) {
                final int dim = parameter.getDimension();

                boolean normalize = (substitutionModel instanceof ComplexSubstitutionModel)
                        && ((ComplexSubstitutionModel) substitutionModel).getNormalization();
                double norm = 1.0;
                double C = 0.0;

                if (normalize) {
                    norm = ((BaseSubstitutionModel) substitutionModel).setupMatrix();
                    double[] qNorm = new double[dim * dim];
                    substitutionModel.getInfinitesimalMatrix(qNorm);
                    for (int a = 0; a < dim; a++) {
                        for (int b = 0; b < dim; b++) {
                            C += gradient[a * dim + b] * qNorm[a * dim + b];
                        }
                    }
                }

                double[] chainedGradient = new double[dim * (dim - 1)];

                int k = 0;
                for (int i = 0; i < dim; ++i) {
                    for (int j = i + 1; j < dim; ++j) {
                        double piJ = parameter.getParameterValue(j);
                        double piI = parameter.getParameterValue(i);
                        chainedGradient[k] = (gradient[i * dim + j] - gradient[i * dim + i] - C * piI) * piJ / norm;
                        ++k;
                    }
                }

                for (int j = 0; j < dim; ++j) {
                    for (int i = j + 1; i < dim; ++i) {
                        double piJ = parameter.getParameterValue(j);
                        double piI = parameter.getParameterValue(i);
                        chainedGradient[k] = (gradient[i * dim + j] - gradient[i * dim + i] - C * piI) * piJ / norm;
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

            @Override
            public Parameter getChainRuleDependent(BastaLikelihood structuredCoalescentLikelihood,
                                                   SubstitutionModel substitutionModel) {
                return substitutionModel.getFrequencyModel().getFrequencyParameter();
            }

        },

        POPULATION_SIZE("populationSize") {
            @Override
            Parameter getParameter(BastaLikelihood structuredCoalescentLikelihood, SubstitutionModel substitutionModel) {
                return structuredCoalescentLikelihood.getPopSizes();
            }

            @Override
            double[] chainRule(double[] gradient, Parameter parameter, SubstitutionModel substitutionModel) {
                final int paramDim = parameter.getDimension();

                if (paramDim == 1 && gradient.length > 1) {
                    // Single shared pop size: sum all per-deme intermediate gradients, then apply chain rule.
                    // By the chain rule: d logL / d N_shared = sum_i (d logL / d N_i).
                    double sum = 0.0;
                    for (double g : gradient) sum += g;
                    double popSize = parameter.getParameterValue(0);
                    return new double[]{sum / -(popSize * popSize)};
                }

                for (int i = 0; i < paramDim; ++i) {
                    double popSize = parameter.getParameterValue(i);
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

            @Override
            public Parameter getChainRuleDependent(BastaLikelihood structuredCoalescentLikelihood,
                                                   SubstitutionModel substitutionModel) {
                return structuredCoalescentLikelihood.getPopSizes();
            }
        };

        WrtParameter(String name) {
            this.name = name;
        }

        abstract Parameter getParameter(BastaLikelihood structuredCoalescentLikelihood, SubstitutionModel substitutionModel);

        abstract double[] chainRule(double[] gradient, Parameter parameter, SubstitutionModel substitutionModel);

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

        public abstract Parameter getChainRuleDependent(BastaLikelihood structuredCoalescentLikelihood,
                                                        SubstitutionModel substitutionModel);
    }
}
