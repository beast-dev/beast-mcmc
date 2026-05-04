/*
 * ProcessGradientSpec.java
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
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

final class ProcessGradientSpec {

    private final List<ContinuousTraitGradientForBranch.ContinuousProcessParameterGradient.DerivationParameter> derivations;

    private ProcessGradientSpec(
            final List<ContinuousTraitGradientForBranch.ContinuousProcessParameterGradient.DerivationParameter> derivations) {
        this.derivations = derivations;
    }

    static ProcessGradientSpec single(
            final ContinuousTraitGradientForBranch.ContinuousProcessParameterGradient.DerivationParameter derivation) {
        return fromList(Collections.singletonList(derivation));
    }

    static ProcessGradientSpec fromList(
            final List<ContinuousTraitGradientForBranch.ContinuousProcessParameterGradient.DerivationParameter> derivations) {
        return new ProcessGradientSpec(
                Collections.unmodifiableList(
                        new ArrayList<ContinuousTraitGradientForBranch.ContinuousProcessParameterGradient.DerivationParameter>(derivations)
                ));
    }

    static ProcessGradientSpec of(
            final ContinuousTraitGradientForBranch.ContinuousProcessParameterGradient.DerivationParameter... derivations) {
        return fromList(Arrays.asList(derivations));
    }

    ContinuousTraitGradientForBranch build(
            final int dim,
            final Tree tree,
            final ContinuousDataLikelihoodDelegate continuousData,
            final String implementation) {
        return ProcessGradientBuildStrategies.build(
                implementation,
                dim,
                tree,
                continuousData,
                derivations);
    }

    List<ContinuousTraitGradientForBranch.ContinuousProcessParameterGradient.DerivationParameter> getDerivations() {
        return derivations;
    }
}
