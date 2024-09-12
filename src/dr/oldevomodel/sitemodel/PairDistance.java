/*
 * PairDistance.java
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

package dr.oldevomodel.sitemodel;

import java.util.List;

/**
 * @author alexei
 */
public class PairDistance implements Comparable<PairDistance> {
    int x, y;
    double distance;
    List<Integer> gaps = null;
    int alignmentLength;

    public PairDistance(int x, int y, double distance, List<Integer> gaps, int alignmentLength) {
        this.x = x;
        this.y = y;
        this.distance = distance;
        this.gaps = gaps;
        this.alignmentLength = alignmentLength;
    }

    public int compareTo(PairDistance p) {
        if (p.distance > distance) return -1;
        if (p.distance < distance) return 1;
        return 0;
    }
}