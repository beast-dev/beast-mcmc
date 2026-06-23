/*
 * CanonicalTransitionMomentProvider.java
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

package dr.evomodel.treedatalikelihood.continuous.canonical;

/**
 * Capability for providers that can expose moment-form transition pieces for
 * an arbitrary effective branch length.
 */
public interface CanonicalTransitionMomentProvider {

    void fillTransitionMatrixFlat(double branchLength, double[] out);

    default void fillTransitionMatrixFlat(final int childNodeIndex,
                                          final double branchLength,
                                          final double[] out) {
        fillTransitionMatrixFlat(branchLength, out);
    }

    void fillTransitionOffset(double branchLength, double[] out);

    default void fillTransitionOffset(final int childNodeIndex,
                                      final double branchLength,
                                      final double[] out) {
        fillTransitionOffset(branchLength, out);
    }

    void fillTransitionCovarianceFlat(double branchLength, double[] out);

    default void fillTransitionCovarianceFlat(final int childNodeIndex,
                                              final double branchLength,
                                              final double[] out) {
        fillTransitionCovarianceFlat(branchLength, out);
    }

    default void fillTransitionMomentsFlat(double branchLength,
                                           double[] transitionMatrixOut,
                                           double[] transitionOffsetOut,
                                           double[] transitionCovarianceOut) {
        fillTransitionMatrixFlat(branchLength, transitionMatrixOut);
        fillTransitionOffset(branchLength, transitionOffsetOut);
        fillTransitionCovarianceFlat(branchLength, transitionCovarianceOut);
    }

    default void fillTransitionMomentsFlat(final int childNodeIndex,
                                           final double branchLength,
                                           final double[] transitionMatrixOut,
                                           final double[] transitionOffsetOut,
                                           final double[] transitionCovarianceOut) {
        fillTransitionMatrixFlat(childNodeIndex, branchLength, transitionMatrixOut);
        fillTransitionOffset(childNodeIndex, branchLength, transitionOffsetOut);
        fillTransitionCovarianceFlat(childNodeIndex, branchLength, transitionCovarianceOut);
    }
}
