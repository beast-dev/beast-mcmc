/*
 * PerBranchCategoryLayoutTest.java
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

package test.dr.evomodel.branchratemodel;

import dr.evomodel.branchratemodel.PerBranchCategoryLayout;
import junit.framework.TestCase;

public class PerBranchCategoryLayoutTest extends TestCase {

    public void testSortsAtomicStatesByReward() {
        // rewards for atomic states 0,1,2,3 are unsorted; sorted order should be 2,0,3,1
        final PerBranchCategoryLayout layout =
                PerBranchCategoryLayout.fromAtomicRewards(new double[]{0.5, 0.9, 0.1, 0.7});

        assertEquals(4, layout.getAtomicStateCount());
        assertEquals(2, layout.getSortedStateIndex(0));
        assertEquals(0, layout.getSortedStateIndex(1));
        assertEquals(3, layout.getSortedStateIndex(2));
        assertEquals(1, layout.getSortedStateIndex(3));
        assertEquals(0.1, layout.getSortedReward(0), 0.0);
        assertEquals(0.9, layout.getSortedReward(3), 0.0);
    }

    public void testInsertionRankBoundaries() {
        final PerBranchCategoryLayout layout =
                PerBranchCategoryLayout.fromAtomicRewards(new double[]{0.2, 0.4, 0.6, 0.8});

        assertEquals("below smallest reward", 0, layout.insertionRank(0.0));
        assertEquals("above largest reward", 4, layout.insertionRank(1.0));
        assertEquals("strictly between two rewards", 2, layout.insertionRank(0.5));
        assertEquals("exactly on a reward value ties to the lower side", 1, layout.insertionRank(0.4));
    }

    public void testBucketCategoryRoundTripAtEveryRank() {
        final double[] rewardByState = {0.5, 0.9, 0.1, 0.7};
        final PerBranchCategoryLayout layout = PerBranchCategoryLayout.fromAtomicRewards(rewardByState);
        final int k = layout.getAtomicStateCount();

        for (int rank = 0; rank <= k; rank++) {
            // continuous category always occupies bucket == rank
            assertEquals(rank, layout.bucketForCategory(rank, PerBranchCategoryLayout.CONTINUOUS_CATEGORY));
            assertEquals(PerBranchCategoryLayout.CONTINUOUS_CATEGORY, layout.categoryForBucket(rank, rank));

            // every bucket decodes to a category that re-encodes to the same bucket
            for (int bucket = 0; bucket <= k; bucket++) {
                final int category = layout.categoryForBucket(rank, bucket);
                assertEquals("bucket " + bucket + " at rank " + rank,
                        bucket, layout.bucketForCategory(rank, category));
            }

            // every atomic category re-decodes to the same category from its bucket
            for (int category = 1; category <= k; category++) {
                final int bucket = layout.bucketForCategory(rank, category);
                assertEquals("category " + category + " at rank " + rank,
                        category, layout.categoryForBucket(rank, bucket));
            }
        }
    }

    public void testInsertionRankZeroPutsContinuousBeforeEveryAtomicState() {
        final PerBranchCategoryLayout layout =
                PerBranchCategoryLayout.fromAtomicRewards(new double[]{0.5, 0.9, 0.1, 0.7});
        final int rank = 0;

        assertEquals(PerBranchCategoryLayout.CONTINUOUS_CATEGORY, layout.categoryForBucket(rank, 0));
        // bucket 1 is now the smallest-reward atomic state (sorted rank 0 -> state 2)
        assertEquals(2 + 1, layout.categoryForBucket(rank, 1));
    }

    public void testInsertionRankAtMaxPutsContinuousAfterEveryAtomicState() {
        final PerBranchCategoryLayout layout =
                PerBranchCategoryLayout.fromAtomicRewards(new double[]{0.5, 0.9, 0.1, 0.7});
        final int k = layout.getAtomicStateCount();

        assertEquals(PerBranchCategoryLayout.CONTINUOUS_CATEGORY, layout.categoryForBucket(k, k));
        // bucket k-1 is the largest-reward atomic state (sorted rank k-1 -> state 1)
        assertEquals(1 + 1, layout.categoryForBucket(k, k - 1));
    }

    public void testSingleAtomicState() {
        final PerBranchCategoryLayout layout = PerBranchCategoryLayout.fromAtomicRewards(new double[]{0.5});

        assertEquals(1, layout.getAtomicStateCount());
        assertEquals(0, layout.insertionRank(0.0));
        assertEquals(1, layout.insertionRank(1.0));
        assertEquals(PerBranchCategoryLayout.CONTINUOUS_CATEGORY, layout.categoryForBucket(0, 0));
        assertEquals(1, layout.categoryForBucket(0, 1));
        assertEquals(1, layout.categoryForBucket(1, 0));
        assertEquals(PerBranchCategoryLayout.CONTINUOUS_CATEGORY, layout.categoryForBucket(1, 1));
    }

    public void testDuplicateRewardValuesStillProduceAValidPermutation() {
        final PerBranchCategoryLayout layout =
                PerBranchCategoryLayout.fromAtomicRewards(new double[]{0.5, 0.5, 0.5});
        final boolean[] seen = new boolean[3];
        for (int rank = 0; rank < 3; rank++) {
            final int state = layout.getSortedStateIndex(rank);
            assertFalse("state " + state + " listed twice", seen[state]);
            seen[state] = true;
        }
        for (boolean s : seen) {
            assertTrue(s);
        }
    }

    public void testRejectsEmptyAtomicStateSet() {
        try {
            PerBranchCategoryLayout.fromAtomicRewards(new double[0]);
            fail("expected IllegalArgumentException");
        } catch (IllegalArgumentException expected) {
            // expected
        }
    }
}
