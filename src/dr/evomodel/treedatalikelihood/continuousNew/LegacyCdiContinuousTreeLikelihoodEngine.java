///*
// * LegacyCdiContinuousTreeLikelihoodEngine.java
// *
// * Copyright © 2002-2024 the BEAST Development Team
// * http://beast.community/about
// *
// * This file is part of BEAST.
// * See the NOTICE file distributed with this work for additional
// * information regarding copyright ownership and licensing.
// *
// * BEAST is free software; you can redistribute it and/or modify
// * it under the terms of the GNU Lesser General Public License as
// * published by the Free Software Foundation; either version 2
// * of the License, or (at your option) any later version.
// *
// *  BEAST is distributed in the hope that it will be useful,
// *  but WITHOUT ANY WARRANTY; without even the implied warranty of
// *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// *  GNU Lesser General Public License for more details.
// *
// * You should have received a copy of the GNU Lesser General Public
// * License along with BEAST; if not, write to the
// * Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
// * Boston, MA  02110-1301  USA
// *
// */
//
//package dr.evomodel.treedatalikelihood.continuousNew;
//
//import dr.evomodel.treedatalikelihood.ProcessOnTreeDelegate.BranchOperation;
//import dr.evomodel.treedatalikelihood.ProcessOnTreeDelegate.NodeOperation;
//
//
//import java.util.List;
//
///**
// * Adapter for the existing CDI-backed continuous tree likelihood path.
// */
//final class LegacyCdiContinuousTreeLikelihoodEngine implements ContinuousTreeLikelihoodEngine {
//
//    private final ContinuousDataLikelihoodDelegateNew delegate;
//
//    LegacyCdiContinuousTreeLikelihoodEngine(final ContinuousDataLikelihoodDelegateNew delegate) {
//        this.delegate = delegate;
//    }
//
//    @Override
//    public double calculateLikelihood(final List<BranchOperation> branchOperations,
//                                      final List<NodeOperation> nodeOperations,
//                                      final int rootNodeNumber) {
//        return delegate.calculateLegacyCdiLikelihood(
//                branchOperations,
//                nodeOperations,
//                rootNodeNumber);
//    }
//
//    @Override
//    public void storeState() {
//        // Legacy CDI state remains owned by ContinuousDataLikelihoodDelegateNew and CDI buffers.
//    }
//
//    @Override
//    public void restoreState() {
//        // Legacy CDI state remains owned by ContinuousDataLikelihoodDelegateNew and CDI buffers.
//    }
//
//    @Override
//    public void acceptState() {
//        // Legacy CDI state remains owned by ContinuousDataLikelihoodDelegateNew and CDI buffers.
//    }
//}
