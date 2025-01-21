/*
 * BiClade.java
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

import java.util.*;

/**
 * @author Andrew Rambaut
 * @version $
 */
class BiClade implements Clade {

    /**
     * Clade for a tip
     * @param index number of the tip
     */
    public BiClade(int index, Taxon taxon) {
        this.index = index;

        count = 0;
        credibility = 1.0;
        size = 1;

        key = index;

        this.taxon = taxon;
    }

    public BiClade(Object key, int size) {
        this.index = -1;
        count = 0;
        credibility = 1.0;
        this.size = size;
        this.key = key;
        this.taxon = null;
    }

    public BiClade(Clade child1, Clade child2) {
        assert child1 instanceof BiClade;
        assert child2 instanceof BiClade;

        count = 0;
        credibility = 0.0;
        BiClade left = (BiClade)child1;
        BiClade right = (BiClade)child2;
        if (right.index < left.index) {
            right = (BiClade)child1;
            left = (BiClade)child2;
        }
        size = left.size + right.size;
        index = left.index;
        addSubClades(left, right);

        key = BiClade.makeKey(left.key, right.key);

        this.taxon = null;
    }

    public void addSubClades(Clade child1, Clade child2) {
        // arrange with the lowest index on the left
        BiClade left = (BiClade)child1;
        BiClade right = (BiClade)child2;
        if (right.index < left.index) {
            right = (BiClade)child1;
            left = (BiClade)child2;
        }
        assert left.size + right.size == size;
        assert left.index == index;
        subClades.add(new Pair<>(left, right));
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
        return taxon;
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

    public Set<Pair<BiClade, BiClade>> getSubClades() {
        return subClades;
    }

    @Override
    public Object getKey() {
        return key;
    }

    public static Object makeKey(Object key1, Object key2) {
        BitSet bits = new BitSet();
        if (key1 instanceof Integer) {
            bits.set((Integer) key1);
        } else {
            assert key1 instanceof BitSet;
            bits.or((BitSet) key1);
        }
        if (key2 instanceof Integer) {
            bits.set((Integer) key2);
        } else {
            assert key2 instanceof BitSet;
            bits.or((BitSet) key2);
        }
        return bits;
    }

//    public boolean equals(Object o) {
//        if (this == o) return true;
//        if (o == null || getClass() != o.getClass()) return false;
//
//        if (((BiClade) o).size != size) return false;
//
//        return !(bits != null ? !Arrays.equals(bits, ((BiClade) o).bits) : ((BiClade) o).bits != null);
//
//    }
//
//    public int hashCode() {
//        return left.hashCode() ^ right.hashCode();
//    }

    public String toString() {
        return "clade " + hashCode();
    }

    private int count;
    private double credibility;
    private final int size;
    private final int index;

    private final Object key;

    private final Taxon taxon;

    private final Set<Pair<BiClade, BiClade>> subClades = new HashSet<>();
    BiClade bestLeft = null;
    BiClade bestRight = null;


    double bestSubTreeCredibility = Double.NaN;
//    double bestSubCladeCredibility = Double.NEGATIVE_INFINITY;

    private List<Object[]> attributeValues = null;
}
