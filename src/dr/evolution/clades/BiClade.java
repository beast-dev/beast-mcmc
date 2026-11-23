/*
 * BiClade.java
 *
 * Copyright © 2002-2026, the BEAST Development Team.
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
import dr.util.Pair;

import java.util.*;
import java.util.concurrent.atomic.AtomicInteger;

/**
 * A clade of a bifurcating tree: a tip, or a set of taxa that can be split into two clades
 * (in one or more ways, if the pairs of sub-clades are being kept).
 *
 * A clade can be counted, given sub-clades, parents and children, and given height and
 * attribute values from any number of threads at once. Its credibility and its best pair of
 * sub-clades are set by whatever summarises the clade system once the counting is done.
 *
 * @author Andrew Rambaut
 * @version $
 */
public class BiClade implements Clade {

    private static final boolean USE_BITSET_CLADE_KEYS = false;

    public static final CladeKeys cladeKeys;

    /**
     * Clade for a tip
     * @param index number of the tip
     */
    public BiClade(int index, Taxon taxon) {
        this.index = index;

        credibility = 1.0;
        size = 1;

        key = cladeKeys.getTaxonKey(index, taxon);

        this.taxon = taxon;
    }

    public BiClade(Object key, int size) {
        this.index = -1;
        credibility = 1.0;
        this.size = size;
        this.key = key;
        this.taxon = null;
    }

    public BiClade(Clade child1, Clade child2) {
        assert child1 instanceof BiClade;
        assert child2 instanceof BiClade;

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

        key = BiClade.getParentKey(left.key, right.key);

        this.taxon = null;
    }

    public synchronized void addSubClades(Clade child1, Clade child2) {
        if (subClades == null) {
            subClades = new HashSet<>();
        }
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

    public synchronized void addParent(BiClade parentClade) {
        if (parentClades == null) {
            parentClades = new HashSet<>();
        }
        parentClades.add(parentClade);
    }

    public synchronized void addChild(BiClade childClade) {
        if (childClades == null) {
            childClades = new HashSet<>();
        }
        childClades.add(childClade);
    }

    public BiClade getMajorityRuleParent() {
        return majorityRuleParent;
    }

    public void setMajorityRuleParent(BiClade majorityRuleParent) {
        this.majorityRuleParent = majorityRuleParent;
    }

    @Override
    public int getCount() {
        return count.get();
    }

    /**
     * Counts one more tree with this clade.
     */
    public void incrementCount() {
        count.incrementAndGet();
    }

    @Override
    public void setCount(int count) {
        this.count.set(count);
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
    public void addAttributeValue(String attributeName, Object value) {
        synchronized (attributeValueMap) {
            attributeValueMap.computeIfAbsent(attributeName, k -> new ArrayList<>()).add(value);
        }
    }

    @Override
    public List<Object> getAttributeValues(String attributeName) {
        return attributeValueMap.get(attributeName);
    }

    @Override
    public void addHeightValue(double height) {
        synchronized (heightValues) {
            heightValues.add(height);
        }
    }

    @Override
    public List<Double> getHeightValues() {
        return heightValues;
    }

//    public void addChildHeightValues(double leftHeight, double rightHeight) {
//        leftHeightValues.add(leftHeight);
//        rightHeightValues.add(rightHeight);
//    }
//
//    public List<Double> getLeftHeightValues() {
//        return leftHeightValues;
//    }
//    public List<Double> getRightHeightValues() {
//        return rightHeightValues;
//    }

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

    public void setBestLeft(BiClade bestLeft) {
        this.bestLeft = bestLeft;
    }

    public void setBestRight(BiClade bestRight) {
        this.bestRight = bestRight;
    }

    public double getBestSubTreeScore() {
        return bestSubTreeScore;
    }

    public void setBestSubTreeScore(double bestSubTreeScore) {
        this.bestSubTreeScore = bestSubTreeScore;
    }

    public synchronized Set<Pair<BiClade, BiClade>> getSubClades() {
        return subClades;
    }

    public synchronized Set<BiClade> getParentClades() {
        return parentClades;
    }

    public synchronized Set<BiClade> getChildClades() {
        return childClades;
    }

    @Override
    public Object getKey() {
        return key;
    }

    public static Object getParentKey(Object key1, Object key2) {
        return cladeKeys.getParentKey(key1, key2);
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (!(o instanceof BiClade)) return false;
        return Objects.equals(key, ((BiClade)o).key);
    }

    @Override
    public int hashCode() {
        return key.hashCode();
    }

    public String toString() {
        return "clade " + key;
    }

    private final AtomicInteger count = new AtomicInteger(0);
    private volatile double credibility;
    private final int size;
    private final int index;

    private final Object key;

    private final Taxon taxon;

    private volatile BiClade majorityRuleParent = null;

    private Set<Pair<BiClade, BiClade>> subClades = null;
    private volatile BiClade bestLeft = null;
    private volatile BiClade bestRight = null;
    private Set<BiClade> parentClades = null;
    private Set<BiClade> childClades = null;


    private volatile double bestSubTreeScore = Double.NaN;

    private final Map<String, List<Object>> attributeValueMap = new HashMap<>();
    private final List<Double> heightValues = new ArrayList<>();

    static {
        if (USE_BITSET_CLADE_KEYS) {
            cladeKeys = BitsetCladeKeys.INSTANCE;
        } else {
            cladeKeys = TaxonHashCladeKeys.INSTANCE;
        }
    }
}
