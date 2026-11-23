/*
 * RunConfigurationTest.java
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

package test.dr.app.beast;

import dr.app.beast.RunConfiguration;
import dr.inference.mcmcmc.MCMCMCOptions;
import junit.framework.TestCase;

import java.util.List;

/**
 * @author Andrew Rambaut
 */
public class RunConfigurationTest extends TestCase {

    private static RunConfiguration configuration(int independentChains, int temperedChains, int chainThreads) {
        double[] temperatures = new double[temperedChains];
        for (int i = 0; i < temperedChains; i++) {
            temperatures[i] = 1.0 / (1.0 + i);
        }
        return new RunConfiguration(independentChains, temperatures, 100,
                MCMCMCOptions.SwapScheme.ORIGINAL_FLAVOR, chainThreads);
    }

    private static String describe(RunConfiguration configuration, int processors) {
        StringBuilder sb = new StringBuilder();
        for (String line : configuration.getThreadingDescription(processors)) {
            sb.append(line).append("\n");
        }
        return sb.toString();
    }

    /**
     * Heated chains count: 16 chains on 4 threads take turns.
     */
    public void testFewerThreadsThanChains() {
        String description = describe(configuration(4, 4, 4), 24);
        assertEquals("Threading:\n" +
                "  4 of the 16 chains compute at a time (--chain_threads 4), on 24 processors.\n" +
                "  The chains take turns, so each computes about 25% of the time. For all of them to compute at once\n" +
                "  use --chain_threads 16, or leave the option out.\n" +
                "  A tempering set advances at the pace of its slowest chain: its chains each run 100 states\n" +
                "  (--mc3_swap) and then wait for one another to swap.\n" +
                "  Each chain has up to 6 threads for its likelihoods while it computes (--thread_scheme auto),\n" +
                "  unless --threads or --beagle_threads is given.\n", description);
    }

    public void testAllChainsAtOnce() {
        String given = describe(configuration(4, 4, 16), 24);
        assertTrue(given, given.contains("  All 16 chains compute at once (--chain_threads 16), on 24 processors.\n"));
        assertFalse(given, given.contains("take turns"));
        assertTrue(given, given.contains("Each chain has one thread for its likelihoods"));

        // the default is every chain, when there are enough processors
        String byDefault = describe(configuration(4, 4, -1), 24);
        assertTrue(byDefault, byDefault.contains("  All 16 chains compute at once, on 24 processors.\n"));
        assertFalse(byDefault, byDefault.contains("take turns"));
    }

    public void testIndependentChainsOnly() {
        String description = describe(configuration(4, 1, 4), 24);
        assertTrue(description, description.contains("  All 4 chains compute at once (--chain_threads 4), on 24 processors.\n"));
        assertFalse(description, description.contains("tempering"));
        assertTrue(description, description.contains("up to 6 threads"));
    }

    /**
     * By default no more chains compute at once than there are processors; there is then
     * nothing to suggest.
     */
    public void testMoreChainsThanProcessors() {
        String description = describe(configuration(8, 4, -1), 8);
        assertTrue(description, description.contains("  8 of the 32 chains compute at a time (one per processor), on 8 processors.\n"));
        assertTrue(description, description.contains("each computes about 25% of the time.\n"));
        assertFalse(description, description.contains("--chain_threads"));
        assertTrue(description, description.contains("Each chain has one thread for its likelihoods"));
    }

    public void testMoreThreadsThanProcessors() {
        String description = describe(configuration(4, 4, 16), 8);
        assertTrue(description, description.contains("  All 16 chains compute at once (--chain_threads 16), on 8 processors.\n"));
        assertTrue(description, description.contains("  More chains compute at once than there are processors, so they share them.\n"));
    }

    public void testOneAtATime() {
        String description = describe(configuration(3, 1, 1), 4);
        assertTrue(description, description.contains("  1 of the 3 chains computes at a time (--chain_threads 1), on 4 processors.\n"));
        assertTrue(description, description.contains("about 33% of the time."));
        assertTrue(description, description.contains("use --chain_threads 3, or leave the option out."));
        assertTrue(description, description.contains("up to 4 threads"));
    }
}
