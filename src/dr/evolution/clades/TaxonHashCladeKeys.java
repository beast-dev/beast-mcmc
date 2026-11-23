/*
 * TaxonHashCladeKeys.java
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

package dr.evolution.clades;

import dr.evolution.util.Taxon;

/**
 * Clade keys that are 64-bit numbers. A taxon's key is a hash of its name and a clade's key
 * is the exclusive-or of the keys of its taxa, so it does not depend on how the clade is
 * built up, on the order or numbering of the taxa, or on anything that varies between runs,
 * threads or chains. There is no state, so this can be used from any thread.
 *
 * Two different clades could in principle have the same key. With 64 bits the chance is
 * negligible for any realistic number of clades.
 *
 * @author Andrew Rambaut
 */
public final class TaxonHashCladeKeys implements CladeKeys {
    public final static CladeKeys INSTANCE = new TaxonHashCladeKeys();

    private TaxonHashCladeKeys() {
    }

    /**
     * @return the key of a taxon: FNV-1a over the characters of its name, then a SplitMix64
     * finaliser so that names differing in a character or two give unrelated keys
     */
    public static long getTaxonKey(String taxonId) {
        long h = 0xcbf29ce484222325L;
        for (int i = 0; i < taxonId.length(); i++) {
            h ^= taxonId.charAt(i);
            h *= 0x100000001b3L;
        }
        h = (h ^ (h >>> 30)) * 0xbf58476d1ce4e5b9L;
        h = (h ^ (h >>> 27)) * 0x94d049bb133111ebL;
        return h ^ (h >>> 31);
    }

    /**
     * @return the key of the clade made of two disjoint clades
     */
    public static long getParentKey(long key1, long key2) {
        return key1 ^ key2;
    }

    @Override
    public Object getTaxonKey(int index, Taxon taxon) {
        return getTaxonKey(taxon.getId());
    }

    @Override
    public Object getParentKey(Object key1, Object key2) {
        return getParentKey(((Long) key1).longValue(), ((Long) key2).longValue());
    }
}
