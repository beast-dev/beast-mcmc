/*
 * IndexSet.java
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

package dr.evomodel.treedatalikelihood.MAT;

/**
 * A set of indices in [0, capacity) that lists its members in the order they were added. Membership
 * is a flag per index, so adding and testing are O(1) and clearing is O(size): for the sets of
 * changed nodes and sites that are filled and emptied once per move.
 */
final class IndexSet {
    private final boolean[] flags;
    private final int[] members;
    private int size = 0;

    IndexSet(int capacity) {
        this.flags = new boolean[capacity];
        this.members = new int[capacity];
    }

    /**
     * @return true if the index was not already a member
     */
    boolean add(int index) {
        if (flags[index]) {
            return false;
        }
        flags[index] = true;
        members[size++] = index;
        return true;
    }

    void addAll(int[] indices, int count) {
        for (int i = 0; i < count; i++) {
            add(indices[i]);
        }
    }

    void addAll(IndexSet other) {
        addAll(other.members, other.size);
    }

    boolean contains(int index) {
        return flags[index];
    }

    int size() {
        return size;
    }

    /**
     * @return the k-th member, in order of addition
     */
    int get(int k) {
        return members[k];
    }

    void clear() {
        for (int i = 0; i < size; i++) {
            flags[members[i]] = false;
        }
        size = 0;
    }
}
