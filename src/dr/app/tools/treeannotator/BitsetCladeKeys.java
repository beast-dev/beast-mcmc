/*
 * BitsetCladeKeys.java
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

/**
 * @author Andrew Rambaut
 * @version $
 */
public class BitsetCladeKeys implements CladeKeys {
    public final static CladeKeys INSTANCE = new BitsetCladeKeys();

    @Override
    public Object getParentKey(Object key1, Object key2) {
        int maxIndex;
        if (key1 instanceof Integer) {
            maxIndex = (Integer) key1;
        } else {
            assert key1 instanceof BitsetKey;
            maxIndex = ((BitsetKey) key1).getMaxIndex();
        }
        if (key2 instanceof Integer) {
            maxIndex = Math.max(maxIndex, (Integer) key2);
        } else {
            assert key2 instanceof BitsetKey;
            maxIndex = Math.max(maxIndex, ((BitsetKey) key2).getMaxIndex());
        }

        BitsetKey key = new BitsetKey(maxIndex);
        if (key1 instanceof Integer) {
            key.set((Integer) key1);
        } else {
            key.setTo((BitsetKey) key1);
        }
        if (key2 instanceof Integer) {
            key.set((Integer) key2);
        } else {
            key.or((BitsetKey) key2);
        }

        return key;
    }

    @Override
    public Object getTaxonKey(int taxon) {
        return taxon;
    }

    private int getMaxIndex(Object key) {
        if (key instanceof Integer) {
            return (Integer) key;
        }
        assert key instanceof BitsetKey;
        return ((BitsetKey) key).getMaxIndex();
    }

    private BitsetKey asBitsetKey(Object key, int maxIndex) {
        if (key instanceof Integer) {
            BitsetKey bitsetKey = new BitsetKey(maxIndex);
            bitsetKey.set((Integer) key);
            return bitsetKey;
        }

        assert key instanceof BitsetKey;
        BitsetKey bitsetKey = new BitsetKey(maxIndex);
        bitsetKey.setTo((BitsetKey) key);
        return bitsetKey;
    }
}
