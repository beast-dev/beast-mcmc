/*
 * DiscreteTraitBranchRateDelegate.java
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

package dr.evomodel.treedatalikelihood.discrete;

import dr.evolution.tree.Tree;
import dr.evomodel.siteratemodel.SiteRateModel;
import dr.evomodel.substmodel.SubstitutionModel;
import dr.evomodel.treedatalikelihood.BeagleDataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.preorder.AbstractBeagleBranchGradientDelegate;

/**
 * @author Xiang Ji
 * @author Marc A. Suchard
 */
public class DiscreteTraitBranchRateDelegate extends AbstractBeagleBranchGradientDelegate {

    static final String GRADIENT_TRAIT_NAME = "BranchRateGradient";
    static final String HESSIAN_TRAIT_NAME = "BranchRateHessian";

    DiscreteTraitBranchRateDelegate(String name, Tree tree, BeagleDataLikelihoodDelegate likelihoodDelegate) {
        super(name, tree, likelihoodDelegate);
    }

    protected void cacheDifferentialMassMatrix(Tree tree, boolean cacheSquaredMatrix) {
        for (int i = 0; i < evolutionaryProcessDelegate.getSubstitutionModelCount(); i++) {
            double[] infinitesimalMatrix = new double[stateCount * stateCount];
            SubstitutionModel substitutionModel = evolutionaryProcessDelegate.getSubstitutionModel(i);
            substitutionModel.getInfinitesimalMatrix(infinitesimalMatrix);

            double[] scaledInfinitesimalMatrix = scaleInfinitesimalMatrixByRates(infinitesimalMatrix, DifferentialChoice.GRADIENT, siteRateModel);
            evolutionaryProcessDelegate.cacheInfinitesimalMatrix(beagle, i, scaledInfinitesimalMatrix);
            if (cacheSquaredMatrix) {
                double[] infinitesimalMatrixSquared = new double[stateCount * stateCount];
                for (int l = 0; l < stateCount; l++) {
                    for (int j = 0; j < stateCount; j++) {
                        double sumOverState = 0.0;
                        for (int k = 0; k < stateCount; k++) {
                            sumOverState += infinitesimalMatrix[l * stateCount + k] * infinitesimalMatrix[k * stateCount + j];
                        }
                        infinitesimalMatrixSquared[l * stateCount + j] = sumOverState;
                    }
                }
                double[] scaledInfinitesimalMatrixSquared = scaleInfinitesimalMatrixByRates(infinitesimalMatrixSquared, DifferentialChoice.HESSIAN, siteRateModel);
                evolutionaryProcessDelegate.cacheInfinitesimalSquaredMatrix(beagle, i, scaledInfinitesimalMatrixSquared);
            }
        }
    }

    static double[] scaleInfinitesimalMatrixByRates(double[] infinitesimalMatrix, DifferentialChoice differentialChoice,
                                                    SiteRateModel siteRateModel) {

        final int matrixSize = infinitesimalMatrix.length;
        double[] scaledInfinitesimalMatrix = new double[matrixSize * siteRateModel.getCategoryCount()];

        for (int i = 0; i < siteRateModel.getCategoryCount(); i++) {
            final double rate = siteRateModel.getRateForCategory(i);

            for (int j = 0; j < matrixSize; j++) {
                scaledInfinitesimalMatrix[matrixSize * i + j] = infinitesimalMatrix[j] * differentialChoice.getRateScale(rate);
            }
        }

        return scaledInfinitesimalMatrix;
    }

    enum DifferentialChoice {
        GRADIENT {
            @Override
            double getRateScale(double rate) {
                return rate;
            }
        },
        HESSIAN {
            @Override
            double getRateScale(double rate) {
                return rate * rate;
            }
        };

        abstract double getRateScale(double rate);
    }

    @SuppressWarnings("unused")
    public static String getName(String name) {
        return GRADIENT_TRAIT_NAME;
    }

    protected String getGradientTraitName() {
        return GRADIENT_TRAIT_NAME;
    }

    protected String getHessianTraitName() {
        return HESSIAN_TRAIT_NAME;
    }

}
