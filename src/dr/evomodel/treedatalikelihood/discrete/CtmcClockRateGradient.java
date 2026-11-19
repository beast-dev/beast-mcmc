/*
 * CtmcClockRateGradient.java
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
import dr.evolution.tree.TreeTrait;
import dr.evomodel.branchratemodel.StrictClockBranchRates;
import dr.evomodel.substmodel.SubstitutionModel;
import dr.evomodel.treedatalikelihood.GradientDataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.inference.hmc.GradientWrtParameterProvider;
import dr.inference.model.Likelihood;
import dr.inference.model.Parameter;
import dr.xml.Reportable;

import java.util.List;

/**
 * Gradient of the log-likelihood with respect to the rate r of a strict clock, from the same cross products
 * (d logL / d Q for each generator Q) that give the gradient with respect to the substitution-model rates.
 * <p>
 * The likelihood depends on r and the generators only through exp(Q r t), so r d logL / d r = sum over the
 * substitution models of sum_ij Q_ij (d logL / d Q)_ij. This holds across rate categories and for an epoch model,
 * whose transition times r does not change. Given a gradient with respect to the rates from the same trait, the
 * clock-rate gradient costs one pass over the S x S cross products per substitution model.
 *
 * @author Marc A. Suchard
 */
public class CtmcClockRateGradient implements GradientWrtParameterProvider, Reportable {

    private final TreeDataLikelihood treeDataLikelihood;
    private final Tree tree;
    private final Parameter rate;
    private final List<SubstitutionModel> substitutionModels;
    private final int stateCount;
    private final SubstitutionModelGradientTraitFactory.Context crossProductTraitContext;
    private final TreeTrait crossProductTrait;

    public CtmcClockRateGradient(String traitName,
                                 TreeDataLikelihood treeDataLikelihood,
                                 GradientDataLikelihoodDelegate likelihoodDelegate,
                                 StrictClockBranchRates branchRates,
                                 AbstractLogAdditiveSubstitutionModelGradient.ApproximationMode mode) {
        this.treeDataLikelihood = treeDataLikelihood;
        this.tree = treeDataLikelihood.getTree();
        this.rate = branchRates.getRateParameter();
        this.substitutionModels = likelihoodDelegate.getBranchModel().getSubstitutionModels();
        this.stateCount = substitutionModels.get(0).getDataType().getStateCount();
        this.crossProductTraitContext = SubstitutionModelGradientTraitFactory.create(traitName,
                treeDataLikelihood,
                likelihoodDelegate,
                stateCount,
                mode == AbstractLogAdditiveSubstitutionModelGradient.ApproximationMode.EXACT_SPECTRAL,
                false,
                mode.getInfo());
        this.crossProductTrait = crossProductTraitContext.getTreeTrait();
    }

    @Override
    public Likelihood getLikelihood() {
        return treeDataLikelihood;
    }

    @Override
    public Parameter getParameter() {
        return rate;
    }

    @Override
    public int getDimension() {
        return 1;
    }

    @Override
    public double[] getGradientLogDensity() {

        crossProductTraitContext.prepareEvaluation();
        final double[] crossProducts = (double[]) crossProductTrait.getTrait(tree, null);

        final int length = stateCount * stateCount;
        final double[] generator = new double[length];

        double total = 0.0;
        for (int m = 0; m < substitutionModels.size(); ++m) {
            substitutionModels.get(m).getInfinitesimalMatrix(generator);
            for (int ij = 0; ij < length; ++ij) {
                total += crossProducts[m * length + ij] * generator[ij];
            }
        }

        return new double[]{total / rate.getParameterValue(0)};
    }

    @Override
    public String getReport() {
        return GradientWrtParameterProvider.getReportAndCheckForError(this, 0.0, Double.POSITIVE_INFINITY, 1E-4);
    }
}
