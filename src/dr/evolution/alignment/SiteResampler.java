/*
 * SiteResampler.java
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

package dr.evolution.alignment;

import java.util.Random;

/**
 * The nonparametric bootstrap's draw: n columns uniformly with replacement from n, as counts per column. Replicates
 * are seeded from a master seed and the replicate number, so any replicate can be drawn again on its own and the
 * draws do not depend on how many replicates there are or in which order they are made.
 *
 * @author Xiang Ji
 * @author Marc A. Suchard
 */
public final class SiteResampler {

    private SiteResampler() {
    }

    /**
     * @return the seed of replicate b: SplitMix64's finalizer of the master seed and b
     */
    public static long replicateSeed(long masterSeed, int replicate) {
        long z = masterSeed + 0x9E3779B97F4A7C15L * (replicate + 1L);
        z = (z ^ (z >>> 30)) * 0xBF58476D1CE4E5B9L;
        z = (z ^ (z >>> 27)) * 0x94D049BB133111EBL;
        return z ^ (z >>> 31);
    }

    /**
     * @return how often each of n columns is drawn when n columns are drawn uniformly with replacement; the counts
     * sum to n
     */
    public static int[] resample(int n, Random random) {
        int[] counts = new int[n];
        for (int i = 0; i < n; i++) {
            counts[random.nextInt(n)]++;
        }
        return counts;
    }
}
