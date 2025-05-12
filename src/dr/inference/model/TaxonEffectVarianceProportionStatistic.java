/*
 * TaxonEffectVarianceProportionStatistic.java
 *
 * Copyright © 2002-2025 the BEAST Development Team
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

package dr.inference.model;

import dr.evolution.tree.Tree;
import dr.evomodel.continuous.MultivariateDiffusionModel;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.evomodel.treedatalikelihood.continuous.TaxonEffectTraitDataModel;
import dr.inference.distribution.DistributionLikelihood;
import dr.inference.distribution.NormalDistributionModel;
import dr.math.matrixAlgebra.Matrix;
import dr.stats.DiscreteStatistics;

/**
 * @author Marc A Suchard
 * @author Philippe Lemey
 */

public class TaxonEffectVarianceProportionStatistic extends AbstractVarianceProportionStatistic
        implements VariableListener, ModelListener {

    private final MultivariateDiffusionModel diffusionModel;
    private final DistributionLikelihood prior;
    private final TreeVarianceSums treeSums;

    private final Parameter effects;

    private Matrix diffusionVariance;
    private Matrix effectsVariance;

    private boolean treeKnown = false;
    private boolean varianceKnown = false;

    private static final boolean USE_POPULATION_VARIANCE = false;

    public TaxonEffectVarianceProportionStatistic(Tree tree, TreeDataLikelihood treeLikelihood,
                                                  TaxonEffectTraitDataModel dataModel,
                                                  MultivariateDiffusionModel diffusionModel,
                                                  DistributionLikelihood prior,
                                                  MatrixRatios ratio) {

        super(tree, treeLikelihood, dataModel, ratio);

        this.diffusionModel = diffusionModel;
        this.treeSums = new TreeVarianceSums(0, 0);
        this.effects = dataModel.getEffects();

        this.diffusionVariance = null;
        this.effectsVariance = null;

        this.prior = prior;

        if (isTreeRandom) {
            ((AbstractModel) tree).addModelListener(this);
        }

        diffusionModel.getPrecisionParameter().addParameterListener(this);
        effects.addParameterListener(this);
    }

    private Matrix getEffectsVariance() {
        // TODO Use variance of prior on effects? [population]

        double variance = DiscreteStatistics.variance(effects.getParameterValues());

        if (USE_POPULATION_VARIANCE &&
                prior != null && prior.getDistribution() instanceof NormalDistributionModel) {
            NormalDistributionModel normal = (NormalDistributionModel) prior.getDistribution();
            variance = normal.variance();
        }
        return new Matrix(new double[][] { new double[] {variance} });
    }

    protected void updateVarianceComponents() {

        double N = tree.getExternalNodeCount();

        double diffusionScale = (treeSums.getDiagonalSum() / N - treeSums.getTotalSum() / (N * N));
        double samplingScale = 1.0; // (N - 1) / N;

        for (int i = 0; i < dimTrait; i++) {

            diffusionComponent.set(i, i, diffusionScale * diffusionVariance.component(i, i));
            samplingComponent.set(i, i, samplingScale * effectsVariance.component(i, i));

            for (int j = i + 1; j < dimTrait; j++) {

                double diffValue = diffusionScale * diffusionVariance.component(i, j);
                double sampValue = samplingScale * effectsVariance.component(i, j);

                diffusionComponent.set(i, j, diffValue);
                samplingComponent.set(i, j, sampValue);

                diffusionComponent.set(j, i, diffValue);
                samplingComponent.set(j, i, sampValue);

            }
        }
    }

    @Override
    protected boolean needToUpdate(int dim) {
        boolean needToUpdate = false;
        if (!treeKnown) {

            updateTreeSums(treeSums);
            treeKnown = true;
            needToUpdate = true;
        }

        if (!varianceKnown) {

            effectsVariance = getEffectsVariance();
//            samplingVariance = taxonEffectTraitDataModel.getSamplingVariance();
            diffusionVariance = new Matrix(diffusionModel.getPrecisionmatrix()).inverse();

            varianceKnown = true;

            needToUpdate = true;
        }
        return needToUpdate;
    }

    @Override
    public void variableChangedEvent(Variable variable, int index, Variable.ChangeType type) {
        assert (variable == effects || variable == diffusionModel.getPrecisionParameter());

        varianceKnown = false;
    }

    @Override
    public void modelChangedEvent(Model model, Object object, int index) {
        assert (model == tree);

        if (!isTreeRandom) throw new IllegalStateException("Attempting to change a fixed tree");

        treeKnown = false;
    }

    @Override
    public void modelRestored(Model model) {
        // Do nothing
    }
}
