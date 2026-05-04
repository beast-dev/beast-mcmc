/*
 * ProcessGradientBuildStrategies.java
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

package dr.evomodelxml.continuous.hmc;

import dr.evolution.tree.Tree;
import dr.evomodel.treedatalikelihood.continuous.ContinuousDataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.continuous.ContinuousTraitGradientForBranch;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

final class ProcessGradientBuildStrategies {

    interface BuildStrategy {
        ContinuousTraitGradientForBranch build(
                int dim,
                Tree tree,
                ContinuousDataLikelihoodDelegate continuousData,
                List<ContinuousTraitGradientForBranch.ContinuousProcessParameterGradient.DerivationParameter> derivations);
    }

    private static final String LEGACY = "legacy";
    private static final String TARGET = "target";

    private static final Map<String, BuildStrategy> STRATEGIES = buildStrategies();

    private ProcessGradientBuildStrategies() {
        // utility
    }

    static ContinuousTraitGradientForBranch build(
            final String strategyName,
            final int dim,
            final Tree tree,
            final ContinuousDataLikelihoodDelegate continuousData,
            final List<ContinuousTraitGradientForBranch.ContinuousProcessParameterGradient.DerivationParameter> derivations) {
        final String normalized = normalizeStrategyName(strategyName);
        final BuildStrategy strategy = STRATEGIES.get(normalized);
        if (strategy == null) {
            throw new IllegalArgumentException("Unknown process-gradient build strategy: " + strategyName);
        }
        return strategy.build(dim, tree, continuousData, derivations);
    }

    private static String normalizeStrategyName(final String strategyName) {
        if (strategyName == null || strategyName.trim().isEmpty()) {
            return LEGACY;
        }
        return strategyName.trim().toLowerCase();
    }

    private static Map<String, BuildStrategy> buildStrategies() {
        final Map<String, BuildStrategy> map = new HashMap<String, BuildStrategy>();

        map.put(LEGACY, new BuildStrategy() {
            @Override
            public ContinuousTraitGradientForBranch build(
                    final int dim,
                    final Tree tree,
                    final ContinuousDataLikelihoodDelegate continuousData,
                    final List<ContinuousTraitGradientForBranch.ContinuousProcessParameterGradient.DerivationParameter> derivations) {
                return new ContinuousTraitGradientForBranch.ContinuousProcessParameterGradient(
                        dim, tree, continuousData, derivations);
            }
        });

        map.put(TARGET, new BuildStrategy() {
            @Override
            public ContinuousTraitGradientForBranch build(
                    final int dim,
                    final Tree tree,
                    final ContinuousDataLikelihoodDelegate continuousData,
                    final List<ContinuousTraitGradientForBranch.ContinuousProcessParameterGradient.DerivationParameter> derivations) {
                final ArrayList<ContinuousTraitGradientForBranch.ProcessGradientTarget> targets =
                        new ArrayList<ContinuousTraitGradientForBranch.ProcessGradientTarget>(derivations.size());
                for (ContinuousTraitGradientForBranch.ContinuousProcessParameterGradient.DerivationParameter derivation : derivations) {
                    targets.add(ContinuousTraitGradientForBranch.targetForDerivationParameter(derivation));
                }
                return new ContinuousTraitGradientForBranch.TargetBasedContinuousProcessParameterGradient(
                        dim, tree, continuousData, targets);
            }
        });

        return Collections.unmodifiableMap(map);
    }
}

