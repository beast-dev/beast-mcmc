/*
 * SubstitutionModelGradientTraitFactory.java
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
import dr.evolution.tree.TreeTraitProvider;
import dr.evomodel.treedatalikelihood.BeagleDataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.DiscreteDataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.GradientDataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.ProcessSimulation;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.evomodel.treedatalikelihood.preorder.ProcessSimulationDelegate;

final class SubstitutionModelGradientTraitFactory {

    private SubstitutionModelGradientTraitFactory() {
    }

    static Context create(String traitName,
                          TreeDataLikelihood treeDataLikelihood,
                          GradientDataLikelihoodDelegate likelihoodDelegate,
                          int stateCount,
                          boolean useExactSpectral,
                          boolean forceAllReal,
                          String defaultInfo) {

        final String treeTraitName = SubstitutionModelCrossProductDelegate.getName(traitName);

        if (treeDataLikelihood.getTreeTrait(treeTraitName) == null) {
            ProcessSimulationDelegate gradientDelegate = createDelegate(traitName,
                    treeDataLikelihood,
                    likelihoodDelegate,
                    stateCount,
                    useExactSpectral,
                    forceAllReal);

            TreeTraitProvider traitProvider = new ProcessSimulation(treeDataLikelihood, gradientDelegate);
            treeDataLikelihood.addTraits(traitProvider.getTreeTraits());
        }

        TreeTrait treeTrait = treeDataLikelihood.getTreeTrait(treeTraitName);
        if (treeTrait == null) {
            throw new IllegalStateException("Could not create tree trait: " + treeTraitName);
        }

        return new Context(treeTrait,
                createPreparer(treeDataLikelihood, likelihoodDelegate),
                getInfo(likelihoodDelegate, defaultInfo));
    }

    private static ProcessSimulationDelegate createDelegate(String traitName,
                                                            TreeDataLikelihood treeDataLikelihood,
                                                            GradientDataLikelihoodDelegate likelihoodDelegate,
                                                            int stateCount,
                                                            boolean useExactSpectral,
                                                            boolean forceAllReal) {

        final Tree tree = treeDataLikelihood.getTree();

        if (likelihoodDelegate instanceof BeagleDataLikelihoodDelegate) {
            BeagleDataLikelihoodDelegate beagleDelegate = (BeagleDataLikelihoodDelegate) likelihoodDelegate;
            if (useExactSpectral) {
                return new SpectralBeagleCrossProductDelegate(traitName,
                        tree,
                        beagleDelegate,
                        treeDataLikelihood.getBranchRateModel(),
                        stateCount);
            }

            return new SubstitutionModelCrossProductDelegate(traitName,
                    tree,
                    beagleDelegate,
                    treeDataLikelihood.getBranchRateModel(),
                    stateCount);
        }

        if (likelihoodDelegate instanceof DiscreteDataLikelihoodDelegate) {
            DiscreteDataLikelihoodDelegate discreteDelegate = (DiscreteDataLikelihoodDelegate) likelihoodDelegate;
            if (discreteDelegate.isSpectralRepresentation()) {
                return new SpectralExactGradientDelegate(traitName,
                        tree,
                        discreteDelegate,
                        stateCount,
                        forceAllReal);
            }

            return new DiscreteSubstitutionModelCrossProductDelegate(traitName,
                    tree,
                    discreteDelegate,
                    stateCount);
        }

        throw new RuntimeException("Other likelihood delegates are currently not supported");
    }

    private static EvaluationPreparer createPreparer(TreeDataLikelihood treeDataLikelihood,
                                                     GradientDataLikelihoodDelegate likelihoodDelegate) {

        if (likelihoodDelegate instanceof DiscreteDataLikelihoodDelegate) {
            final DiscreteDataLikelihoodDelegate discreteDelegate =
                    (DiscreteDataLikelihoodDelegate) likelihoodDelegate;

            return new EvaluationPreparer() {
                @Override
                public void prepare() {
                    if (!discreteDelegate.hasKnownBranchLengths()) {
                        discreteDelegate.updatePostOrdersFromTreeDataLikelihood(treeDataLikelihood);
                    }
                }
            };
        }

        return EvaluationPreparer.NONE;
    }

    private static String getInfo(GradientDataLikelihoodDelegate likelihoodDelegate,
                                  String defaultInfo) {
        if (likelihoodDelegate instanceof DiscreteDataLikelihoodDelegate
                && ((DiscreteDataLikelihoodDelegate) likelihoodDelegate).isSpectralRepresentation()) {
            return "an exact spectral Frechet derivative";
        }

        return defaultInfo;
    }

    interface EvaluationPreparer {
        EvaluationPreparer NONE = new EvaluationPreparer() {
            @Override
            public void prepare() {
            }
        };

        void prepare();
    }

    static final class Context {
        private final TreeTrait treeTrait;
        private final EvaluationPreparer preparer;
        private final String info;

        private Context(TreeTrait treeTrait,
                        EvaluationPreparer preparer,
                        String info) {
            this.treeTrait = treeTrait;
            this.preparer = preparer;
            this.info = info;
        }

        TreeTrait getTreeTrait() {
            return treeTrait;
        }

        void prepareEvaluation() {
            preparer.prepare();
        }

        String getInfo() {
            return info;
        }
    }
}
