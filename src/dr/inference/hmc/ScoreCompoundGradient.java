/*
 * ScoreCompoundGradient.java
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

package dr.inference.hmc;

import java.util.ArrayList;
import java.util.List;

/**
 * A {@link CompoundGradient} that also exposes the additive local score rows of its components,
 * concatenated in component order, so an empirical Fisher information can be assembled in one
 * analytical pass.  Components that do not provide score rows contribute zero columns; if no
 * component provides any, {@link #getScoreContributionsLogDensity()} returns {@code null}.
 *
 * <p>This is a separate subclass, rather than a change to {@link CompoundGradient}, so the widely
 * used compound gradient is unaffected.</p>
 */
public class ScoreCompoundGradient extends CompoundGradient implements ScoreContributionsWrtParameterProvider {

    public ScoreCompoundGradient(List<GradientWrtParameterProvider> derivativeList) {
        super(derivativeList);
    }

    @Override
    public double[][] getScoreContributionsLogDensity() {
        final List<double[][]> parts = new ArrayList<double[][]>(derivativeList.size());
        int rowCount = -1;
        for (final GradientWrtParameterProvider provider : derivativeList) {
            if (!(provider instanceof ScoreContributionsWrtParameterProvider)) {
                parts.add(null);
                continue;
            }
            final double[][] contribution =
                    ((ScoreContributionsWrtParameterProvider) provider).getScoreContributionsLogDensity();
            if (contribution == null) {
                parts.add(null);
                continue;
            }
            if (rowCount < 0) {
                rowCount = contribution.length;
            } else if (contribution.length != rowCount) {
                throw new IllegalStateException("Compound score contributions have unequal row counts.");
            }
            parts.add(contribution);
        }
        if (rowCount < 0) {
            return null;
        }
        final double[][] result = new double[Math.max(0, rowCount)][dimension];
        int offset = 0;
        for (int partIndex = 0; partIndex < parts.size(); ++partIndex) {
            final int partDimension = derivativeList.get(partIndex).getDimension();
            final double[][] part = parts.get(partIndex);
            if (part == null) {
                offset += partDimension;
                continue;
            }
            for (int row = 0; row < result.length; ++row) {
                if (part[row] == null || part[row].length != partDimension) {
                    throw new IllegalStateException("Compound score contribution dimension mismatch.");
                }
                System.arraycopy(part[row], 0, result[row], offset, partDimension);
            }
            offset += partDimension;
        }
        return result;
    }
}
