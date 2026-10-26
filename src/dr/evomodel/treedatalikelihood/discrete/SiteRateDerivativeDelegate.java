/*
 * SiteRateDerivativeDelegate.java
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

import dr.evolution.tree.Tree;
import dr.evolution.tree.BranchRates;
import dr.evomodel.siteratemodel.DifferentiableSiteRateModel;
import dr.evomodel.substmodel.SubstitutionModel;
import dr.evomodel.treedatalikelihood.BeagleDataLikelihoodDelegate;
import dr.inference.model.Parameter;

/**
 * Per-branch derivatives of log L with respect to the length of each branch, where the rate categories move together
 * as d r_c / d theta for a site-rate parameter theta: BEAGLE edge differentials with (d r_c / d theta) Q~ cached in
 * place of the branch-rate delegate's r_c Q~. Summing length_e times the result over branches gives the part of
 * d log L / d theta that flows through the category rates.
 * <p>
 * The infinitesimal-matrix buffers are shared with the branch-rate and node-height delegates of the same likelihood,
 * which cache r_c Q~ and may treat them as still valid, so this delegate restores r_c Q~ after each use.
 *
 * @author Xiang Ji
 * @author Marc A. Suchard
 */
final class SiteRateDerivativeDelegate extends DiscreteTraitBranchRateDelegate {

    private static final String GRADIENT_TRAIT_NAME = "SiteRateGradient";
    private static final String HESSIAN_TRAIT_NAME = "SiteRateHessian";

    private final String traitName;
    private final DifferentiableSiteRateModel differentiableSiteRateModel;
    private final Parameter parameter;

    SiteRateDerivativeDelegate(String name, Tree tree, BeagleDataLikelihoodDelegate likelihoodDelegate,
                               BranchRates branchRateModel, DifferentiableSiteRateModel siteRateModel,
                               Parameter parameter) {
        super(name, tree, likelihoodDelegate, branchRateModel);
        this.traitName = name;
        this.differentiableSiteRateModel = siteRateModel;
        this.parameter = parameter;
    }

    @Override
    protected void cacheDifferentialMassMatrix(Tree tree, boolean cacheSquaredMatrix) {
        final int categoryCount = differentiableSiteRateModel.getCategoryCount();
        final double[] dRates = new double[categoryCount];
        differentiableSiteRateModel.getCategoryDifferentials(parameter, dRates, new double[categoryCount]);

        final int matrixSize = stateCount * stateCount;
        for (int i = 0; i < evolutionaryProcessDelegate.getSubstitutionModelCount(); i++) {
            final double[] infinitesimalMatrix = new double[matrixSize];
            final SubstitutionModel substitutionModel = evolutionaryProcessDelegate.getSubstitutionModel(i);
            substitutionModel.getInfinitesimalMatrix(infinitesimalMatrix);

            final double[] scaled = new double[matrixSize * categoryCount];
            for (int c = 0; c < categoryCount; c++) {
                for (int j = 0; j < matrixSize; j++) {
                    scaled[matrixSize * c + j] = infinitesimalMatrix[j] * dRates[c];
                }
            }
            evolutionaryProcessDelegate.cacheInfinitesimalMatrix(beagle, i, scaled);
        }
    }

    @Override
    protected void getNodeDerivatives(Tree tree, double[] first, double[] second) {
        if (second != null) {
            throw new UnsupportedOperationException("Second derivatives wrt site-rate parameters are not implemented");
        }
        substitutionProcessKnown = false; // another delegate may have written the shared buffers
        super.getNodeDerivatives(tree, first, null);
        super.cacheDifferentialMassMatrix(tree, false); // put back r_c Q~ for the branch-rate and node-height delegates
        substitutionProcessKnown = false;
    }

    static String gradientTraitName(String name) {
        return GRADIENT_TRAIT_NAME + ":" + name;
    }

    @Override
    protected String getGradientTraitName() {
        return gradientTraitName(traitName);
    }

    @Override
    protected String getHessianTraitName() {
        return HESSIAN_TRAIT_NAME + ":" + traitName;
    }
}
