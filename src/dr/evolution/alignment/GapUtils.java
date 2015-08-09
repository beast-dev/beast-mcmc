/*
 * GapUtils.java
 *
 * Copyright (c) 2002-2015 Alexei Drummond, Andrew Rambaut and Marc Suchard
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
 */

package dr.evolution.alignment;

import java.util.List;

/**
 * @author Alexei Drummond
 *
 * @version $Id: GapUtils.java,v 1.3 2005/04/20 21:26:19 rambaut Exp $
 */
public class GapUtils {

    /**
     * Adds Integer objects to the given list representing the sizes of the gaps found in the given alignment.
     * Note -- it only handles pairwise alignments
     * @param alignment a pairwise alignment
     * @param gaps a list in which the gap sizes will be added
     */
    public static void getGapSizes(Alignment alignment, List gaps) {

        int stateCount = alignment.getDataType().getStateCount();
        int gapStart = Integer.MAX_VALUE;
        for (int i = 0; i < alignment.getSiteCount(); i++) {
            int x = alignment.getState(0,i);
            int y = alignment.getState(1,i);

            if (y < stateCount && x < stateCount) {
                // no gap
                if (gapStart < i) {
                    gaps.add(new Integer(i-gapStart));
                    gapStart = Integer.MAX_VALUE;
                }
            } else {
                // gap
                if (gapStart > i) gapStart = i;
            }
        }
    }
}
