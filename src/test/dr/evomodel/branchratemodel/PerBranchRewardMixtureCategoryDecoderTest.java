/*
 * PerBranchRewardMixtureCategoryDecoderTest.java
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

import dr.evomodel.branchratemodel.PerBranchRewardMixtureCategoryDecoder;
import dr.evomodel.branchratemodel.RewardRates;
import dr.inference.model.Parameter;
import junit.framework.TestCase;

public class PerBranchRewardMixtureCategoryDecoderTest extends TestCase {

    // 4 atomic states with rewards (unsorted by state index): 0->0.7, 1->0.1, 2->0.9, 3->0.3
    // sorted ascending: state 1 (0.1), state 3 (0.3), state 0 (0.7), state 2 (0.9)
    private static RewardRates fourStateRewardRates() {
        return new RewardRates(
                new Parameter.Default("values", new double[]{0.7, 0.1, 0.9, 0.3}),
                null,
                new Parameter.Default("varying", new double[]{0.7, 0.1, 0.9, 0.3}),
                new Parameter.Default("stateIndices", new double[]{0.0, 1.0, 2.0, 3.0}));
    }

    private static PerBranchRewardMixtureCategoryDecoder decoderForCts(final double... ctsValues) {
        final int dim = ctsValues.length;
        final Parameter categoryParameter = new Parameter.Default("category", new double[dim]);
        for (int i = 0; i < dim; i++) {
            categoryParameter.setParameterValue(i, 0.5); // arbitrary, inside category-0's [0,1] bucket
        }
        final Parameter cuts = new Parameter.Default("cuts", new double[]{0, 1, 2, 3, 4, 5});
        final Parameter cts = new Parameter.Default("cts", ctsValues);
        return new PerBranchRewardMixtureCategoryDecoder(categoryParameter, cuts, cts, fourStateRewardRates(), 4, dim);
    }

    public void testContinuousBranchDecodesToCategoryZero() {
        // cts=0.05 is below every atomic reward, so insertion rank is 0 and bucket 0
        // (where the raw coordinate 0.5 decodes to) is the continuous category.
        final PerBranchRewardMixtureCategoryDecoder decoder = decoderForCts(0.05);
        assertEquals(0, decoder.getCategoryForParameterIndex(0));
        assertFalse(decoder.isAtomic(0));
    }

    public void testAtomicBucketMapsToSortedNeighbourNotRawStateIndex() {
        // cts = 0.05, below every atomic reward -> insertion rank 0 -> continuous at bucket 0,
        // bucket 1 is the smallest-reward atomic state (state 1, reward 0.1), NOT state 0.
        final PerBranchRewardMixtureCategoryDecoder decoder = decoderForCts(0.05);
        final Parameter category = decoder.getCategoryParameter();
        category.setParameterValue(0, 1.5); // bucket 1
        assertEquals(1 + 1, decoder.getCategoryForParameterIndex(0)); // category = state 1 + 1
        assertTrue(decoder.isAtomic(0));
        assertEquals(1, decoder.getAtomicState(0));
    }

    public void testDifferentBranchesGetDifferentNeighboursFromTheSameRawCoordinate() {
        // two branches, cts values on opposite ends of the atomic-reward range
        final PerBranchRewardMixtureCategoryDecoder decoder = decoderForCts(0.05, 0.95);
        final Parameter category = decoder.getCategoryParameter();
        category.setParameterValue(0, 1.5); // bucket 1
        category.setParameterValue(1, 1.5); // same raw bucket, different branch

        // branch 0 (cts=0.05, rank 0): bucket 1 -> smallest-reward atomic state (state 1)
        assertEquals(1, decoder.getAtomicState(0));
        // branch 1 (cts=0.95, rank 4): bucket 1 -> second-smallest-reward atomic state (state 3)
        assertEquals(3, decoder.getAtomicState(1));
    }

    public void testGetCategoryForValueMatchesParameterIndexDecoding() {
        final PerBranchRewardMixtureCategoryDecoder decoder = decoderForCts(0.05, 0.95);
        assertEquals(decoder.getCategoryForParameterIndex(0), decoder.getCategoryForValue(0, 0.5));
        assertEquals(decoder.getCategoryForParameterIndex(1), decoder.getCategoryForValue(1, 0.5));
    }

    public void testRefreshEmbeddingPicksUpChangedCtsValue() {
        final PerBranchRewardMixtureCategoryDecoder decoder = decoderForCts(0.05);
        final Parameter category = decoder.getCategoryParameter();
        category.setParameterValue(0, 1.5); // bucket 1
        assertEquals(1, decoder.getAtomicState(0)); // rank 0: bucket 1 -> state 1 (smallest reward)

        decoder.getCutParameter(); // no-op, just exercising the accessor
        // move this branch's cts above every atomic reward and refresh
        // (cts is the "cts" Parameter passed to decoderForCts; grab it back out via reflection-free
        // trick: rebuild via a fresh decoder instead, since the fixture doesn't expose it directly)
        final Parameter categoryParameter = new Parameter.Default("category2", new double[]{1.5});
        final Parameter cuts = new Parameter.Default("cuts2", new double[]{0, 1, 2, 3, 4, 5});
        final Parameter cts = new Parameter.Default("cts2", new double[]{0.05});
        final PerBranchRewardMixtureCategoryDecoder movable =
                new PerBranchRewardMixtureCategoryDecoder(categoryParameter, cuts, cts, fourStateRewardRates(), 4, 1);
        assertEquals(1, movable.getAtomicState(0));

        cts.setParameterValue(0, 0.95);
        movable.refreshEmbedding();
        assertEquals(3, movable.getAtomicState(0)); // rank 4: bucket 1 -> second-smallest-reward state
    }

    public void testLowerAndUpperCutAreBranchSpecific() {
        final PerBranchRewardMixtureCategoryDecoder decoder = decoderForCts(0.05, 0.95);
        // category 0 (continuous) sits at bucket = insertion rank, which differs per branch
        assertEquals(0.0, decoder.getLowerCut(0, 0), 0.0);
        assertEquals(1.0, decoder.getUpperCut(0, 0), 0.0);
        assertEquals(4.0, decoder.getLowerCut(1, 0), 0.0);
        assertEquals(5.0, decoder.getUpperCut(1, 0), 0.0);
    }

    public void testRejectsMismatchedCtsDimension() {
        final Parameter categoryParameter = new Parameter.Default("category", new double[]{0.5, 0.5});
        final Parameter cuts = new Parameter.Default("cuts", new double[]{0, 1, 2, 3, 4, 5});
        final Parameter cts = new Parameter.Default("cts", new double[]{0.5}); // wrong dimension
        try {
            new PerBranchRewardMixtureCategoryDecoder(categoryParameter, cuts, cts, fourStateRewardRates(), 4, 2);
            fail("expected IllegalArgumentException");
        } catch (IllegalArgumentException expected) {
            // expected
        }
    }
}
