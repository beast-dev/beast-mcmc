/*
 * RuntimeOverridesTest.java
 *
 * Copyright © 2002-2026 the BEAST Development Team
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

package test.dr.inference.mcmc;

import dr.inference.mcmc.RuntimeOverrides;
import org.junit.After;
import org.junit.Test;

import static org.junit.Assert.*;

/**
 * The chain length and logging frequency given in a file, and the command line options that override them.
 */
public class RuntimeOverridesTest {

    @After
    public void clearProperties() {
        System.clearProperty(RuntimeOverrides.CHAIN_LENGTH_PROPERTY);
        System.clearProperty(RuntimeOverrides.LOG_EVERY_PROPERTY);
    }

    @Test
    public void theFileIsFollowedWhenNothingIsOverridden() {
        assertEquals(1000000L, RuntimeOverrides.getChainLength(1000000L));
        assertEquals(10000, RuntimeOverrides.getLogEvery(10000));
    }

    @Test
    public void anOverriddenChainLengthIsUsed() {
        System.setProperty(RuntimeOverrides.CHAIN_LENGTH_PROPERTY, "5000");
        assertEquals(5000L, RuntimeOverrides.getChainLength(1000000L));
    }

    @Test
    public void anOverriddenLoggingFrequencyIsUsed() {
        System.setProperty(RuntimeOverrides.LOG_EVERY_PROPERTY, "250");
        assertEquals(250, RuntimeOverrides.getLogEvery(10000));
    }

    @Test
    public void aLongerChainThanTheFileAsksForIsAllowed() {
        System.setProperty(RuntimeOverrides.CHAIN_LENGTH_PROPERTY, "20000000");
        assertEquals(20000000L, RuntimeOverrides.getChainLength(1000L));
    }

    @Test
    public void aLoggerWithNoFixedFrequencyIsLeftAlone() {
        System.setProperty(RuntimeOverrides.LOG_EVERY_PROPERTY, "250");
        // a tree logger that was not given a logEvery, or one that logs at every state
        assertEquals(0, RuntimeOverrides.getLogEvery(0));
        assertEquals(-1, RuntimeOverrides.getLogEvery(-1));
    }

    @Test
    public void aValueThatIsNotAWholeNumberIsIgnored() {
        System.setProperty(RuntimeOverrides.CHAIN_LENGTH_PROPERTY, "banana");
        System.setProperty(RuntimeOverrides.LOG_EVERY_PROPERTY, "1e6");

        assertEquals("the file's value should stand", 1000000L, RuntimeOverrides.getChainLength(1000000L));
        assertEquals(10000, RuntimeOverrides.getLogEvery(10000));
    }

    @Test
    public void aValueThatIsNotPositiveIsIgnored() {
        System.setProperty(RuntimeOverrides.CHAIN_LENGTH_PROPERTY, "0");
        System.setProperty(RuntimeOverrides.LOG_EVERY_PROPERTY, "-100");

        assertEquals(1000000L, RuntimeOverrides.getChainLength(1000000L));
        assertEquals(10000, RuntimeOverrides.getLogEvery(10000));
    }

    @Test
    public void aLoggingFrequencyBeyondAnIntIsCappedRatherThanWrapped() {
        System.setProperty(RuntimeOverrides.LOG_EVERY_PROPERTY, Long.toString(4L * Integer.MAX_VALUE));
        assertEquals(Integer.MAX_VALUE, RuntimeOverrides.getLogEvery(10000));
    }
}
