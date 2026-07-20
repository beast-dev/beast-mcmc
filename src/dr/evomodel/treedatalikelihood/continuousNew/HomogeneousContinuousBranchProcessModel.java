/*
 * HomogeneousContinuousBranchProcessModel.java
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

import dr.evomodel.continuous.canonical.ContinuousProcessModel;

/**
 * One process model applied to every non-root branch.
 */
public final class HomogeneousContinuousBranchProcessModel implements ContinuousBranchProcessModel {

    private final ContinuousProcessModel processModel;
    private final ContinuousBranchTimeModel branchTimeModel;

    public HomogeneousContinuousBranchProcessModel(final ContinuousProcessModel processModel,
                                                   final ContinuousBranchTimeModel branchTimeModel) {
        this.processModel = processModel;
        this.branchTimeModel = branchTimeModel;
    }

    @Override
    public int getProcessModelCount() {
        return 1;
    }

    @Override
    public ContinuousProcessModel getProcessModel(final int processIndex) {
        if (processIndex != 0) {
            throw new IllegalArgumentException("Homogeneous continuous process index must be 0: " + processIndex);
        }
        return processModel;
    }

    @Override
    public int getSegmentCount(final int childNodeIndex) {
        return branchTimeModel.getSegmentCount(childNodeIndex);
    }

    @Override
    public int getSegmentProcessIndex(final int childNodeIndex, final int segmentIndex) {
        branchTimeModel.getSegmentProcessIndex(childNodeIndex, segmentIndex);
        return 0;
    }

    @Override
    public double getSegmentEffectiveTime(final int childNodeIndex, final int segmentIndex) {
        return branchTimeModel.getSegmentEffectiveTime(childNodeIndex, segmentIndex);
    }
}
