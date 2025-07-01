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

    void addParent(BiClade parentClade) {
        if (parentClades == null) {
            parentClades = new HashSet<>();
        }
        parentClades.add(parentClade);
    }

    void addChild(BiClade childClade) {
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
        synchronized (attributeValues) {
            attributeValues.add(values);
        }
    }

    @Override
    public List<Object[]> getAttributeValues() {
        return attributeValues;
    }

    public void addHeightValue(double height) {
        synchronized (heightValues) {
            heightValues.add(height);
        }
    }

    public List<Double> getHeightValues() {
        return heightValues;
    }
    public void addChildHeightValues(double leftHeight, double rightHeight) {
        leftHeightValues.add(leftHeight);
        rightHeightValues.add(rightHeight);
    }

    public List<Double> getLeftHeightValues() {
        return leftHeightValues;
    }
    public List<Double> getRightHeightValues() {
        return rightHeightValues;
    }

    public void setMeanHeight(double meanHeight) {
        this.meanHeight = meanHeight;
    }

    public void setMedianHeight(double medianHeight) {
        this.medianHeight = medianHeight;
    }

    public void setHeightRange(Double[] range) {
        this.heightRange = range;
    }

    public void setHeightHPD(Double[] HPDs) {
        this.heightHPDs = HPDs;
    }

    public double getMeanHeight() {
        return meanHeight;
    }

    public double getMedianHeight() {
        return medianHeight;
    }

    public Double[] getHeightRange() {
        return heightRange;
    }

    public Double[] getHeightHPDs() {
        return heightHPDs;
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

    public Set<BiClade> getParentClades() {
        return parentClades;
    }

    public Set<BiClade> getChildClades() {
        return childClades;
    }

    @Override
    public Object getKey() {
        return key;
    }

    public static Object makeKey(Object key1, Object key2) {
        int maxIndex;
        if (key1 instanceof Integer) {
            maxIndex = (Integer) key1;
        } else {
            assert key1 instanceof CladeKey;
            maxIndex = ((CladeKey) key1).getMaxIndex();
        }
        if (key2 instanceof Integer) {
            maxIndex = Math.max(maxIndex, (Integer) key2);
        } else {
            assert key2 instanceof CladeKey;
            maxIndex = Math.max(maxIndex, ((CladeKey) key2).getMaxIndex());
        }

        CladeKey key = new CladeKey(maxIndex);
        if (key1 instanceof Integer) {
            key.set((Integer) key1);
        } else {
            key.setTo((CladeKey) key1);
        }
        if (key2 instanceof Integer) {
            key.set((Integer) key2);
        } else {
            key.or((CladeKey) key2);
        }

        return key;
    }
    public static Object makeKey(Object... keys) {
        int maxIndex = 0;
        for (Object key : keys) {
            maxIndex = Math.max(maxIndex, key instanceof Integer ? (Integer) key : ((CladeKey) key).getMaxIndex());
        }
        CladeKey bits = new CladeKey(maxIndex);
        for (Object key : keys) {
            if (key instanceof Integer) {
                bits.set((Integer) key);
            } else {
                assert key instanceof CladeKey;
                bits.or((CladeKey) key);
            }
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

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (!(o instanceof BiClade)) return false;
        return Objects.equals(key, ((BiClade)o).key);
    }

    @Override
    public int hashCode() {
        return key.hashCode();
//        return Objects.hash(key);
    }

    public String toString() {
        return "clade " + key;
    }

    int count;
    double credibility;
    final int size;
    final int index;

    final Object key;

    private final Taxon taxon;

    public BiClade majorityRuleParent = null;

    private Set<Pair<BiClade, BiClade>> subClades = null;
    BiClade bestLeft = null;
    BiClade bestRight = null;
    private Set<BiClade> parentClades = null;
    private Set<BiClade> childClades = null;


    double bestSubTreeScore = Double.NaN;

    private final List<Object[]> attributeValues = new ArrayList<>();
    private final List<Double> heightValues = new ArrayList<>();
    private final List<Double> leftHeightValues = new ArrayList<>();
    private final List<Double> rightHeightValues = new ArrayList<>();

    private double meanHeight;
    private double medianHeight;
    private Double[] heightRange;
    private Double[] heightHPDs ;

}
