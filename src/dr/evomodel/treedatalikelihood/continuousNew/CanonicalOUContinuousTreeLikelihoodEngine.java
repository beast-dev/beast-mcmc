/*
 * CanonicalOUContinuousTreeLikelihoodEngine.java
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

package dr.evomodel.treedatalikelihood.continuousNew;

import dr.evomodel.treedatalikelihood.ProcessOnTreeDelegate.BranchOperation;
import dr.evomodel.treedatalikelihood.ProcessOnTreeDelegate.NodeOperation;
import dr.evomodel.treedatalikelihood.continuousNew.canonical.adapter.CanonicalOUIntegrator;

import java.util.List;

/**
 * Adapter for the existing canonical OU tree likelihood backend.
 */
final class CanonicalOUContinuousTreeLikelihoodEngine implements ContinuousTreeLikelihoodEngine {

    private final ContinuousDataLikelihoodDelegateNew delegate;
    private final CanonicalOUIntegrator canonicalOUIntegrator;

    CanonicalOUContinuousTreeLikelihoodEngine(final ContinuousDataLikelihoodDelegateNew delegate,
                                              final CanonicalOUIntegrator canonicalOUIntegrator) {
        this.delegate = delegate;
        this.canonicalOUIntegrator = canonicalOUIntegrator;
    }

    @Override
    public double calculateLikelihood(final List<BranchOperation> branchOperations,
                                      final List<NodeOperation> nodeOperations,
                                      final int rootNodeNumber) {
        return delegate.calculateCanonicalOULikelihood();
    }

    @Override
    public void storeState() {
        canonicalOUIntegrator.storeState();
    }

    @Override
    public void restoreState() {
        canonicalOUIntegrator.restoreState();
    }

    @Override
    public void acceptState() {
        canonicalOUIntegrator.acceptState();
    }
}
