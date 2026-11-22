/*
 * ArgumentsSectionTest.java
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

package test.dr.app.util;

import dr.app.util.Arguments;
import org.junit.Test;

import java.io.ByteArrayOutputStream;
import java.io.PrintStream;

import static org.junit.Assert.*;

/**
 * Headings that group the options in the usage. A heading is not an option: it is listed in the usage but
 * cannot be given on the command line.
 */
public class ArgumentsSectionTest {

    private Arguments arguments() {
        return new Arguments(new Arguments.Option[]{
                new Arguments.Section("First group"),
                new Arguments.Option("alpha", "a", "the first option"),
                new Arguments.Section("Second group"),
                new Arguments.IntegerOption("beta", "b", "the second option"),
        });
    }

    private String usageOf(Arguments arguments) {
        PrintStream out = System.out;
        ByteArrayOutputStream captured = new ByteArrayOutputStream();
        try {
            System.setOut(new PrintStream(captured));
            arguments.printUsage("tool", "<file>");
        } finally {
            System.setOut(out);
        }
        return captured.toString();
    }

    @Test
    public void theHeadingsAppearInTheUsage() {
        String usage = usageOf(arguments());

        assertTrue(usage, usage.contains("First group:"));
        assertTrue(usage, usage.contains("Second group:"));
        assertTrue(usage, usage.contains("--alpha | -a: the first option"));
        assertTrue(usage, usage.contains("--beta | -b: the second option"));
    }

    @Test
    public void aHeadingIsNotListedAmongTheOptionsToGive() {
        String firstLine = usageOf(arguments()).split("\n")[0];

        assertTrue(firstLine, firstLine.contains("[--alpha|-a]"));
        assertFalse("the summary line should list options only", firstLine.contains("First group"));
    }

    @Test
    public void theOptionsStillParse() throws Exception {
        Arguments arguments = arguments();
        arguments.parseArguments(new String[]{"-alpha", "-beta", "7", "data.xml"});

        assertTrue(arguments.hasOption("alpha"));
        assertEquals(7, arguments.getIntegerOption("beta"));
        assertEquals(1, arguments.getLeftoverArguments().length);
        assertEquals("data.xml", arguments.getLeftoverArguments()[0]);
    }

    @Test
    public void aHeadingCannotBeGivenAsAnOption() {
        Arguments arguments = arguments();
        try {
            arguments.parseArguments(new String[]{"-First group", "data.xml"});
        } catch (Arguments.ArgumentException expected) {
            // an unknown option is the right answer; either way it must not be treated as a heading
        }
        assertFalse(arguments.hasOption("First group"));
    }
}
