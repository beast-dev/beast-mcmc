/*
 * FastClade.java
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

package dr.app.tools.treeannotator;

import dr.evolution.util.Taxon;
import dr.util.Pair;

import java.math.BigInteger;
import java.security.MessageDigest;
import java.security.NoSuchAlgorithmException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Set;

/**
 * @author Andrew Rambaut
 * @version $
 */
class FastClade implements Clade {
    public FastClade(int index) {
        this.index = index;

        count = 0;
        credibility = 1.0;
        size = 1;
        bits = null;
        hash = index;
        key = index;
    }

    public FastClade(Clade leftClade, Clade rightClade, int tipCount) {
        assert leftClade instanceof FastClade;
        assert rightClade instanceof FastClade;

        FastClade left = (FastClade)leftClade;
        FastClade right = (FastClade)rightClade;

        count = 0;
        credibility = 0.0;
        size = left.size + right.size;
        index = -1;

        bits = new byte[(tipCount / 8) + 1];
        if (left.bits == null) {
            int byteIndex = left.index / 8;
            int bitMask = 1 << (left.index % 8);
            bits[byteIndex] = (byte) bitMask;
        } else {
            System.arraycopy(left.bits, 0, bits, 0, left.bits.length);
        }

        if (right.bits == null) {
            int byteIndex = right.index / 8;
            int bitMask = 1 << (right.index % 8);
            bits[byteIndex] |= (byte) bitMask;
        } else {
            for (int i = 0; i < bits.length; i++) {
                bits[i] |= right.bits[i];
            }
        }

//            FastClade sc1, sc2;
//            if (subClade1.index < subClade2.index) {
//                index = subClade1.index;
//                sc1 = subClade1;
//                sc2 = subClade2;
//            } else {
//                index = subClade2.index;
//                sc1 = subClade2;
//                sc2 = subClade1;
//            }
//
//            int maxIndex = Math.max(sc1.maxIndex(), sc2.maxIndex());
//            bits = new byte[maxIndex - index + 1];
//
//            if (sc1.bits == null) {
//                bits[0] = 0b10000000;
//            } else {
//                System.arraycopy(sc1.bits, 0, bits, 0, sc1.bits.length);
//            }
//
//            if (sc2.bits == null) {
//                bits[sc2.index - index] = 1;
//            } else {
//                for (int i = 0; i < sc2.bits.length; i++) {
//                    bits[i + sc2.index - index] |= sc2.bits[i];
//                }
//            }

        MessageDigest digest = null;
        try {
            digest = MessageDigest.getInstance("SHA-1");
            byte[] encodedhash = digest.digest(bits);
            key = new BigInteger(1, encodedhash);
            //hash = Arrays.hashCode(encodedhash);
            hash = 0;
        } catch (NoSuchAlgorithmException e) {
            throw new RuntimeException(e);
        }

    }

    private int maxIndex() {
        if (bits != null) {
            return (index + bits.length - 1);
        } else {
            return index;
        }
    }

    private int byteIndex(int index) {
        return index >> 8;
    }

    private int bitIndex(int index) {
        return index | 8;
    }

    @Override
    public int getCount() {
        return count;
    }

    @Override
    public void setCount(int count) {
        this.count = count;
    }

    @Override
    public double getCredibility() {
        return credibility;
    }

    @Override
    public void setCredibility(double credibility) {
        this.credibility = credibility;
    }

    @Override
    public void addAttributeValues(Object[] values) {
        if (attributeValues == null) {
            attributeValues = new ArrayList<>();
        }
        attributeValues.add(values);
    }

    @Override
    public List<Object[]> getAttributeValues() {
        return attributeValues;
    }

    @Override
    public int getSize() {
        return size;
    }

    @Override
    public Taxon getTaxon() {
        return null;
    }

    @Override
    public int getIndex() {
        return index;
    }

    @Override
    public Clade getBestLeft() {
        return bestLeft;
    }

    @Override
    public Clade getBestRight() {
        return bestRight;
    }

    public Integer getHash() {
        return hash;
    }

    @Override
    public Object getKey() {
        return key;
    }

    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        if (((FastClade) o).size != size) return false;

        return !(bits != null ? !Arrays.equals(bits, ((FastClade) o).bits) : ((FastClade) o).bits != null);

    }

    public int hashCode() {
        return hash;
    }

    public String toString() {
        return "clade " + hashCode();
    }

    int count;
    double credibility;
    final int size;
    final byte[] bits;
    final int index;

    final Object key;
    final int hash;
    List<Object[]> attributeValues = null;
    Set<Pair<FastClade, FastClade>> subClades = null;
    FastClade bestLeft = null;
    FastClade bestRight = null;
    double bestSubTreeCredibility;
}
