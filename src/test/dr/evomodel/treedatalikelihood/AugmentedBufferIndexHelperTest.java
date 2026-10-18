/*
 * AugmentedBufferIndexHelperTest.java
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

package test.dr.evomodel.treedatalikelihood;

import dr.evomodel.treedatalikelihood.AugmentedBufferIndexHelper;
import dr.evomodel.treedatalikelihood.BufferIndexHelper;
import junit.framework.TestCase;

import java.util.HashSet;
import java.util.Random;
import java.util.Set;

/**
 * Checks AugmentedBufferIndexHelper in the three layouts that use it (partials, partials with pre-order slots, and
 * matrices with cached infinitesimal matrices) against BufferIndexHelper, the layout it replaces.
 *
 * @author Marc A Suchard
 */
public class AugmentedBufferIndexHelperTest extends TestCase {

    // 6 tips and 3 transition times (4 epochs)
    private static final int TIPS = 6;
    private static final int NODES = 2 * TIPS - 1;
    private static final int CAPACITY = TIPS * 3;
    private static final int CACHED = 4 * 4; // K = 4E: Q and Q^2, double-buffered, for each epoch

    private static final int[] ALLOCATED = {0, 1, 7, CAPACITY};

    public AugmentedBufferIndexHelperTest(String name) {
        super(name);
    }

    /**
     * One use of the helper. The fixed buffers after the copies of the original nodes (pre-order partials of the
     * original nodes, or cached matrices) are [2 * NODES - min, blockBase). Today, the same buffers were
     * BufferIndexHelper(NODES + CAPACITY, min) followed by everything else.
     */
    private static final class Layout {
        final String name;
        final int min;
        final int blockBase;
        final int stride;
        final int todayCount; // buffers allocated today, at CAPACITY

        Layout(String name, int min, int blockBase, int stride, int todayCount) {
            this.name = name;
            this.min = min;
            this.blockBase = blockBase;
            this.stride = stride;
            this.todayCount = todayCount;
        }

        AugmentedBufferIndexHelper create() {
            return new AugmentedBufferIndexHelper(NODES, min, CAPACITY, blockBase, stride);
        }

        int fixedStart() {
            return 2 * NODES - min;
        }
    }

    private static Layout[] layouts() {
        final int todayPartials = new BufferIndexHelper(NODES + CAPACITY, TIPS).getBufferCount();
        final int todayMatrices = new BufferIndexHelper(NODES + CAPACITY, 0).getBufferCount();
        return new Layout[]{
                new Layout("partials", TIPS, 2 * NODES - TIPS, 2, todayPartials),
                new Layout("partials with pre-order", TIPS, 3 * NODES - TIPS, 3,
                        todayPartials + NODES + CAPACITY),
                new Layout("matrices", 0, 2 * NODES + CACHED, 2, todayMatrices + CACHED),
        };
    }

    /**
     * Every buffer of A augmented nodes (both copies, other slots) and of the original nodes is distinct, and
     * together with the fixed buffers they are exactly 0 ... getBufferCount() - 1. At CAPACITY that is today's set.
     */
    public void testIndicesDistinctAndBelowCount() {
        for (Layout layout : layouts()) {
            for (int a : ALLOCATED) {
                AugmentedBufferIndexHelper helper = layout.create();
                helper.setAllocated(a);
                final int count = helper.getBufferCount();
                assertEquals(layout.blockBase + layout.stride * a, count);
                assertEquals(count, helper.getBufferCount(a));

                Set<Integer> indices = new HashSet<Integer>();
                for (int id = 0; id < NODES + a; ++id) {
                    final int copy0 = helper.getOffsetIndex(id);
                    add(indices, copy0, count, layout.name);
                    if (id >= layout.min) {
                        helper.flipOffset(id);
                        add(indices, helper.getOffsetIndex(id), count, layout.name);
                    }
                    if (id >= NODES) {
                        assertEquals(helper.getSlotIndex(id, 0), copy0);
                        assertEquals(helper.getSlotIndex(id, 1), helper.getOffsetIndex(id));
                        for (int slot = 2; slot < layout.stride; ++slot) {
                            add(indices, helper.getSlotIndex(id, slot), count, layout.name);
                        }
                    }
                }
                for (int index = layout.fixedStart(); index < layout.blockBase; ++index) {
                    add(indices, index, count, layout.name);
                }
                assertEquals(layout.name + ", A = " + a, count, indices.size());

                if (a == CAPACITY) {
                    assertEquals(layout.name, layout.todayCount, count);
                }
            }
        }
    }

    private static void add(Set<Integer> indices, int index, int count, String name) {
        assertTrue(name + ": index " + index, index >= 0 && index < count);
        assertTrue(name + ": index " + index + " used twice", indices.add(index));
    }

    /**
     * Allocating more augmented nodes changes no index, whatever the flips, stores and restores before it.
     */
    public void testNoIndexChangesWhenGrowing() {
        Random random = new Random(1);
        for (Layout layout : layouts()) {
            AugmentedBufferIndexHelper helper = layout.create();
            helper.setAllocated(ALLOCATED[0]);

            for (int k = 1; k < ALLOCATED.length; ++k) {
                final int before = ALLOCATED[k - 1];
                for (int n = 0; n < 50; ++n) {
                    apply(helper, randomOperation(random, layout.min, NODES + before), Integer.MAX_VALUE);
                }

                int[] indices = new int[NODES + before];
                boolean[] safe = new boolean[NODES + before];
                for (int id = 0; id < NODES + before; ++id) {
                    indices[id] = helper.getOffsetIndex(id);
                    safe[id] = id >= layout.min && helper.isSafeUpdate(id);
                }

                helper.setAllocated(ALLOCATED[k]);

                for (int id = 0; id < NODES + before; ++id) {
                    assertEquals(layout.name + ", id " + id, indices[id], helper.getOffsetIndex(id));
                    assertEquals(safe[id], id >= layout.min && helper.isSafeUpdate(id));
                }
            }
        }
    }

    /**
     * Under random flips, stores and restores, the original nodes have the indices of BufferIndexHelper(NODES, min),
     * and every id follows today's BufferIndexHelper(NODES + CAPACITY, min): the same copy and the same isSafeUpdate.
     */
    public void testFlipStoreRestoreMatchToday() {
        Random random = new Random(2);
        for (Layout layout : layouts()) {
            AugmentedBufferIndexHelper helper = layout.create();
            BufferIndexHelper fixed = new BufferIndexHelper(NODES, layout.min);
            BufferIndexHelper today = new BufferIndexHelper(NODES + CAPACITY, layout.min);

            for (int step = 0; step < 2000; ++step) {
                final int operation = randomOperation(random, layout.min, NODES + CAPACITY);
                apply(helper, operation, Integer.MAX_VALUE);
                apply(fixed, operation, NODES);
                apply(today, operation, Integer.MAX_VALUE);

                for (int id = 0; id < NODES + CAPACITY; ++id) {
                    final int index = helper.getOffsetIndex(id);
                    if (id < NODES) {
                        assertEquals(layout.name + ", id " + id, fixed.getOffsetIndex(id), index);
                    }
                    if (id < layout.min) {
                        assertEquals(id, today.getOffsetIndex(id));
                        continue;
                    }

                    final int copy = today.getOffsetIndex(id) == id ? 0 : 1;
                    final int expected = (id < NODES) ? id + copy * (NODES - layout.min) :
                            helper.getSlotIndex(id, copy);
                    assertEquals(layout.name + ", step " + step + ", id " + id, expected, index);
                    assertEquals(today.isSafeUpdate(id), helper.isSafeUpdate(id));
                }
            }
        }
    }

    private static final int STORE = -1;
    private static final int RESTORE = -2;

    /**
     * @return a flip of an id in [min, end) (the id itself), a store or a restore
     */
    private static int randomOperation(Random random, int min, int end) {
        final int choice = random.nextInt(10);
        return (choice < 7) ? min + random.nextInt(end - min) : ((choice < 9) ? STORE : RESTORE);
    }

    /**
     * Applies an operation; a helper that covers only ids below end ignores flips of the others.
     */
    private static void apply(BufferIndexHelper helper, int operation, int end) {
        if (operation == STORE) {
            helper.storeState();
        } else if (operation == RESTORE) {
            helper.restoreState();
        } else if (operation < end) {
            helper.flipOffset(operation);
        }
    }
}
