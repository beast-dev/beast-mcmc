/*
 * CanonicalPreparedBranchSnapshot.java
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

package dr.evomodel.treedatalikelihood.continuousNew.canonical;

import dr.evomodel.continuous.ou.canonical.CanonicalPreparedBranchHandle;
import dr.evomodel.treedatalikelihood.continuousNew.canonical.message.CanonicalGaussianTransition;

/**
 * Prepared, branch-local canonical transition data for one child node.
 */
public final class CanonicalPreparedBranchSnapshot {

    private int childNodeIndex;
    private int processIndex;
    private double effectiveBranchLength;
    private CanonicalGaussianTransition transition;
    private CanonicalPreparedBranchHandle preparedBranchHandle;

    public CanonicalPreparedBranchSnapshot(final int childNodeIndex,
                                           final double effectiveBranchLength,
                                           final CanonicalGaussianTransition transition,
                                           final CanonicalPreparedBranchHandle preparedBranchHandle) {
        this(childNodeIndex, 0, effectiveBranchLength, transition, preparedBranchHandle);
    }

    public CanonicalPreparedBranchSnapshot(final int childNodeIndex,
                                           final int processIndex,
                                           final double effectiveBranchLength,
                                           final CanonicalGaussianTransition transition,
                                           final CanonicalPreparedBranchHandle preparedBranchHandle) {
        this.childNodeIndex = childNodeIndex;
        this.processIndex = processIndex;
        this.effectiveBranchLength = effectiveBranchLength;
        this.transition = transition;
        this.preparedBranchHandle = preparedBranchHandle;
    }

    public void reset(final int childNodeIndex,
                      final double effectiveBranchLength,
                      final CanonicalGaussianTransition transition,
                      final CanonicalPreparedBranchHandle preparedBranchHandle) {
        reset(childNodeIndex, 0, effectiveBranchLength, transition, preparedBranchHandle);
    }

    public void reset(final int childNodeIndex,
                      final int processIndex,
                      final double effectiveBranchLength,
                      final CanonicalGaussianTransition transition,
                      final CanonicalPreparedBranchHandle preparedBranchHandle) {
        this.childNodeIndex = childNodeIndex;
        this.processIndex = processIndex;
        this.effectiveBranchLength = effectiveBranchLength;
        this.transition = transition;
        this.preparedBranchHandle = preparedBranchHandle;
    }

    public int getChildNodeIndex() {
        return childNodeIndex;
    }

    public int getProcessIndex() {
        return processIndex;
    }

    public double getEffectiveBranchLength() {
        return effectiveBranchLength;
    }

    public CanonicalGaussianTransition getTransition() {
        return transition;
    }

    public CanonicalPreparedBranchHandle getPreparedBranchHandle() {
        return preparedBranchHandle;
    }
}
