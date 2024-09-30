/*
 * AbstractGlmSubstitutionModelGradient.java
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

import dr.evomodel.substmodel.GlmSubstitutionModel;
import dr.evomodel.treedatalikelihood.BeagleDataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.inference.distribution.GeneralizedLinearModel;
import dr.inference.loggers.LogColumn;
import dr.inference.loggers.Loggable;
import dr.inference.model.CompoundParameter;
import dr.inference.model.Parameter;
import dr.util.Author;
import dr.util.Citation;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

/**
 * @author Marc A. Suchard
 */

public abstract class AbstractGlmSubstitutionModelGradient extends AbstractLogAdditiveSubstitutionModelGradient {

    protected final GeneralizedLinearModel glm;
    private final ParameterMap parameterMap;

    public AbstractGlmSubstitutionModelGradient(String traitName,
                                                TreeDataLikelihood treeDataLikelihood,
                                                BeagleDataLikelihoodDelegate likelihoodDelegate,
                                                GlmSubstitutionModel substitutionModel,
                                                ApproximationMode mode) {

        super(traitName, treeDataLikelihood, likelihoodDelegate, substitutionModel, mode);
        this.glm = substitutionModel.getGeneralizedLinearModel();
        this.parameterMap = makeParameterMap(glm);
    }

    @SuppressWarnings("unused")
    String getType() { return "fixed"; }

    ParameterMap makeParameterMap(GeneralizedLinearModel glm) {

        final List<Integer> whichBlock = new ArrayList<>();
        final List<Integer> whichIndex = new ArrayList<>();

        CompoundParameter cp = new CompoundParameter("fixedEffects");
        boolean multi = glm.getNumberOfFixedEffects() > 1;

        for (int i = 0; i < glm.getNumberOfFixedEffects(); ++i) {
            Parameter p = glm.getFixedEffect(i);
            if (multi) {
                cp.addParameter(p);
            }
            for (int j = 0; j < p.getDimension(); ++j) {
                whichBlock.add(i);
                whichIndex.add(j);
            }

            if (glm.getFixedEffectIndicator(i) != null) {
                throw new IllegalArgumentException("GLM fixed effects gradients do not currently work with indicator variables");
            }
        }

        final Parameter whichParameter = multi ? cp : glm.getFixedEffect(0);
        return new ParameterMap() {

            @Override
            public double[] getCovariateColumn(int i) {
                return glm.getDesignMatrix(whichBlock.get(i)).getColumnValues(whichIndex.get(i));
            }

            @Override
            public Parameter getParameter() { return whichParameter; }
        };
    }

    @Override
    public Parameter getParameter() {
        return parameterMap.getParameter();
    }

    protected double preProcessNormalization(double[] differentials, double[] generator,
                                             boolean normalize) {
        return 0.0;
    }

    double processSingleGradientDimension(int i,
                                          double[] differentials, double[] generator, double[] pi,
                                          boolean normalize, double normalizationConstant) {

        double[] covariate = parameterMap.getCovariateColumn(i);
        return calculateCovariateDifferential(generator, differentials, covariate, pi, normalize);
    }

    private double calculateCovariateDifferential(double[] generator, double[] crossProduct,
                                                  double[] covariate, double[] pi,
                                                  boolean doNormalization) {

        double normalization = 0.0;
        double total = 0.0;

        int k = 0;
        for (int i = 0; i < stateCount; ++i) {
            for (int j = i + 1; j < stateCount; ++j) {

                double xij = covariate[k++];
                double element = xij * generator[index(i,j)];

                if (element != 0.0) {
                    total += crossProduct[index(i, j)] * element;
                    total -= crossProduct[index(i, i)] * element;

                    total += correction(i, j, crossProduct) * element;
                    total -= correction(i, i, crossProduct) * element;

                    normalization += element * pi[i];
                }
            }
        }

        for (int j = 0; j < stateCount; ++j) {
            for (int i = j + 1; i < stateCount; ++i) {

                double xij = covariate[k++];
                double element = xij * generator[index(i,j)];

                if (element != 0.0) {
                    total += crossProduct[index(i, j)] * element;
                    total -= crossProduct[index(i, i)] * element;

                    total += correction(i, j, crossProduct) * element;
                    total -= correction(i, j, crossProduct) * element;

                    normalization += element * pi[i];
                }
            }
        }

        if (doNormalization) {
            for (int i = 0; i < stateCount; ++i) {
                for (int j = 0; j < stateCount; ++j) {
                    total -= crossProduct[index(i,j)] * generator[index(i,j)] * normalization;
                }
            }
        }

        return total;
    }

    int index(int i, int j) {
        return i * stateCount + j;
    }

    @Override
    public LogColumn[] getColumns() {
        return Loggable.getColumnsFromReport(this, "gradient report");
    }

    @Override
    public Citation.Category getCategory() {
        return Citation.Category.FRAMEWORK;
    }

    @Override
    public String getDescription() {
        return "Using linear-time differential calculations for all substitution generator elements";
    }

    @Override
    public List<Citation> getCitations() {
        return Collections.singletonList(CITATION);
    }

    private static final Citation CITATION = new Citation(
            new Author[]{
                    new Author( "AF", "Magee"),
                    new Author( "AJ", "Holbrook"),
                    new Author( "JE", "Pekar"),
                    new Author( "IW", "Caviedes-Solis"),
                    new Author( "FA", "Matsen"),
                    new Author( "G", "Baele"),
                    new Author( "JO", "Wertheim"),
                    new Author( "X", "Ji"),
                    new Author("P", "Lemey"),
                    new Author("MA", "Suchard"),
            },
            "Random-effects substitution models for phylogenetics via scalable gradient approximations",
            "",
            Citation.Status.IN_PREPARATION);

    interface ParameterMap {
        double[] getCovariateColumn(int i);
        Parameter getParameter();
    }
}
