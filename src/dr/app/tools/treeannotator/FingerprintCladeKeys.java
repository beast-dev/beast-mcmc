/*
 * FingerprintCladeKeys.java
 *
 * Copyright © 2026, the BEAST Development Team.
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
 */

package dr.app.tools.treeannotator;

import dr.math.MathUtils;

import java.util.*;

/**
 * Implements clade keys using the random tip numbers and XOR compositing method of Patrick Varilly
 * described in https://doi.org/10.1101/2025.03.25.645253
 *
 * @author Andrew Rambaut
 * @version $
 */
public final class FingerprintCladeKeys implements CladeKeys {
    public final static CladeKeys INSTANCE = new FingerprintCladeKeys();

    private FingerprintCladeKeys() {}

    @Override
    public Object getParentKey(Object key1, Object key2) {
        long fingerprint = ((Long)key1) ^ ((Long)key2);
        return fingerprint;
    }

    @Override
    public Object getTaxonKey(int taxon) {
        long fingerprint = TAXON_FINGERPRINTS.computeIfAbsent(taxon, integer -> {
            long fp = MathUtils.nextLong();
            // abundance of caution - check each fingerprint is unique
            while (FINGERPRINTS.contains(fp)) {
                assert true : "Taxon fingerprint collision"; // if exceptions are on then flag and stop
                fp = MathUtils.nextLong();
            }
            FINGERPRINTS.add(fp);
            return fp;
        });
        TAXON_FINGERPRINTS.put(taxon, fingerprint);
        return fingerprint;
    }

    private final static Map<Integer, Long> TAXON_FINGERPRINTS = new HashMap<>();
    private final static Set<Long> FINGERPRINTS = new HashSet<>();
}
