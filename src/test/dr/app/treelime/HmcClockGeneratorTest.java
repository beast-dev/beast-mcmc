/*
 * HmcClockGeneratorTest.java
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

package test.dr.app.treelime;

import dr.app.treelime.options.TreeLimeOptions;
import dr.app.treelime.options.PartitionClockModel;
import dr.app.treelime.types.ClockType;
import junit.framework.TestCase;

/**
 * The HMC relaxed clock XML that TreeLime writes, with one clock and with two unlinked (prefixed) clocks, must parse
 * with every idref resolved and each clock's priors in the prior block.
 */
public class HmcClockGeneratorTest extends TestCase {

    public HmcClockGeneratorTest(String name) {
        super(name);
    }

    public void testOneClock() throws Exception {
        check(false);
    }

    public void testTwoUnlinkedClocks() throws Exception {
        check(true);
    }

    private void check(boolean unlinkClocks) throws Exception {
        TreeLimeOptions options = TreeLimeHarness.importData(TreeLimeHarness.TWO_PARTITIONS);
        if (options == null) {
            return;
        }
        if (unlinkClocks) {
            TreeLimeHarness.unlinkClockModels(options);
            assertEquals(2, options.getPartitionClockModels().size());
        }
        for (PartitionClockModel clock : options.getPartitionClockModels()) {
            clock.setClockType(ClockType.HMC_CLOCK);
        }
        String xml = TreeLimeHarness.generateAndParse(options);
        for (PartitionClockModel clock : options.getPartitionClockModels()) {
            String prefix = clock.getPrefix();
            for (String prior : new String[]{"ratesPrior", "scalePrior"}) {
                assertTrue("prior block refers to " + prefix + prior, xml.contains("idref=\"" + prefix + prior + "\""));
            }
        }
    }
}
