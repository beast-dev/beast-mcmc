/*
 * MixedEffectsClockGeneratorTest.java
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

package test.dr.app.beauti;

import dr.app.beauti.options.BeautiOptions;
import dr.app.beauti.options.PartitionClockModel;
import dr.app.beauti.types.ClockType;
import dr.evolution.util.Taxa;
import junit.framework.TestCase;

/**
 * The mixed-effects clock XML that BEAUti writes must parse (plan 3, E7): one fixed effect and one category per
 * monophyletic taxon set, consecutively numbered, every id prefixed so that several clock models can coexist, and
 * each clock's priors in the prior block.
 */
public class MixedEffectsClockGeneratorTest extends TestCase {

    public MixedEffectsClockGeneratorTest(String name) {
        super(name);
    }

    public void testOneClockWithNonMonophyleticTaxonSet() throws Exception {
        check(false);
    }

    public void testTwoUnlinkedClocks() throws Exception {
        check(true);
    }

    private void check(boolean unlinkClocks) throws Exception {
        BeautiOptions options = BeautiGeneratorHarness.importData(BeautiGeneratorHarness.TWO_PARTITIONS);
        if (options == null) {
            return;
        }

        // a non-monophyletic set first, so the numbering of the monophyletic ones is exercised
        int[][] sets = {{0, 7}, {1, 2}, {3, 4}};
        boolean[] mono = {false, true, true};
        for (int i = 0; i < sets.length; i++) {
            Taxa taxa = new Taxa("set" + (i + 1));
            for (int index : sets[i]) {
                taxa.addTaxon(options.taxonList.getTaxon(index));
            }
            options.taxonSets.add(taxa);
            options.taxonSetsMono.put(taxa, mono[i]);
            options.taxonSetsIncludeStem.put(taxa, i == 2);
            options.taxonSetsHeights.put(taxa, 0.0);
            options.taxonSetsTreeModel.put(taxa, options.getPartitionTreeModels().get(0));
        }

        if (unlinkClocks) {
            BeautiGeneratorHarness.unlinkClockModels(options);
            assertEquals(2, options.getPartitionClockModels().size());
        }
        for (PartitionClockModel clock : options.getPartitionClockModels()) {
            clock.setClockType(ClockType.MIXED_EFFECTS_CLOCK);
        }

        String xml = BeautiGeneratorHarness.generateAndParse(options);
        assertTrue(xml.contains("usePreOrder=\"true\""));
        for (PartitionClockModel clock : options.getPartitionClockModels()) {
            String prefix = clock.getPrefix();
            // the rate prior is left out when BEAUti fixes the rate (no tip dates in this data set)
            for (String prior : new String[]{"scalePrior", "interceptPrior",
                    "fixedEffectsLikelihood1", "fixedEffectsLikelihood2"}) {
                assertTrue("prior block refers to " + prefix + prior,
                        xml.contains("idref=\"" + prefix + prior + "\""));
            }
        }
    }
}
