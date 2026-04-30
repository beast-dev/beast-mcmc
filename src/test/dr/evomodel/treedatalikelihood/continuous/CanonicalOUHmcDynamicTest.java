/*
 * CanonicalOUHmcDynamicTest.java
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

package test.dr.evomodel.treedatalikelihood.continuous;

import junit.framework.TestCase;

import java.io.File;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.TimeUnit;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class CanonicalOUHmcDynamicTest extends TestCase {

    private static final String REPO_HMC_XML =
            "src/test/resources/dr/evomodel/treedatalikelihood/continuous/canonical_ou_hmc_seed666.xml";
    private static final long REGRESSION_CHAIN_LENGTH = 100L;
    private static final double HMC_FINAL_STATE_TOLERANCE = 1.5e-2;
    private static final Pattern SCREEN_LOG_ROW = Pattern.compile(
            "(?m)^(100)\\s+([-+0-9.Ee]+)\\s+([-+0-9.Ee]+)\\s+([-+0-9.Ee]+)\\s+" +
                    "([-+0-9.Ee]+)\\s+([-+0-9.Ee]+)\\s+.*$");

    public CanonicalOUHmcDynamicTest(final String name) {
        super(name);
    }

    public void testSeed666CanonicalOrthogonalHmcFinalState() throws Exception {
        final File xml = new File(REPO_HMC_XML).getCanonicalFile();
        assertTrue("Missing repo-owned canonical OU HMC regression fixture: " + xml, xml.isFile());

        final File shortenedXml = createShortenedHmcFixture(xml);
        try {
            final HmcScreenLogRow finalRow = runBeastAndReadFinalScreenLogRow(shortenedXml);

            assertEquals("state", REGRESSION_CHAIN_LENGTH, finalRow.state);
            assertEquals("joint", -937.3183, finalRow.joint, HMC_FINAL_STATE_TOLERANCE);
            assertEquals("prior", -1125.1585, finalRow.prior, HMC_FINAL_STATE_TOLERANCE);
            assertEquals("likelihood", 187.8402, finalRow.likelihood, HMC_FINAL_STATE_TOLERANCE);
            assertEquals("age(root)", 32.7070, finalRow.rootAge, 5.0e-5);
            assertEquals("ucgd.mean", 1.00000, finalRow.ucgdMean, 5.0e-6);
        } finally {
            if (!shortenedXml.delete()) {
                shortenedXml.deleteOnExit();
            }
        }
    }

    private HmcScreenLogRow runBeastAndReadFinalScreenLogRow(final File xml) throws Exception {
        final File output = File.createTempFile("canonical-ou-hmc-seed666", ".out");
        try {
            final Process process = new ProcessBuilder(buildBeastCommand(xml))
                    .directory(new File(".").getCanonicalFile())
                    .redirectErrorStream(true)
                    .redirectOutput(output)
                    .start();

            if (!process.waitFor(45, TimeUnit.SECONDS)) {
                process.destroyForcibly();
                fail("Timed out running BEAST HMC regression fixture");
            }

            final String screenLog = readFile(output);
            assertEquals("BEAST process exit code\n" + tail(screenLog), 0, process.exitValue());
            return parseFinalRow(screenLog);
        } finally {
            if (!output.delete()) {
                output.deleteOnExit();
            }
        }
    }

    private List<String> buildBeastCommand(final File xml) {
        final List<String> command = new ArrayList<String>();
        command.add(new File(System.getProperty("java.home"), "bin/java").getPath());
        command.add("-cp");
        command.add(System.getProperty("java.class.path"));
        command.add("dr.app.beast.BeastMain");
        command.add("-seed");
        command.add("666");
        command.add("-overwrite");
        command.add("-citations_off");
        command.add(xml.getPath());
        return command;
    }

    private File createShortenedHmcFixture(final File sourceXml) throws IOException {
        final String xml = readFile(sourceXml);
        final String shortened = xml.replace(
                "<mcmc id=\"mcmc\" chainLength=\"1000\" autoOptimize=\"true\">",
                "<mcmc id=\"mcmc\" chainLength=\"" + REGRESSION_CHAIN_LENGTH + "\" autoOptimize=\"true\">");
        assertFalse("Expected to shorten the HMC regression fixture", xml.equals(shortened));

        final File target = File.createTempFile("canonical-ou-hmc-seed666-chain100", ".xml");
        Files.write(target.toPath(), shortened.getBytes(StandardCharsets.UTF_8));
        return target;
    }

    private HmcScreenLogRow parseFinalRow(final String screenLog) {
        final Matcher matcher = SCREEN_LOG_ROW.matcher(screenLog);
        HmcScreenLogRow row = null;
        while (matcher.find()) {
            row = new HmcScreenLogRow(
                    Long.parseLong(matcher.group(1)),
                    Double.parseDouble(matcher.group(2)),
                    Double.parseDouble(matcher.group(3)),
                    Double.parseDouble(matcher.group(4)),
                    Double.parseDouble(matcher.group(5)),
                    Double.parseDouble(matcher.group(6)));
        }
        assertNotNull("Could not find final state 100 in BEAST screen log\n" + tail(screenLog), row);
        return row;
    }

    private static String readFile(final File file) throws IOException {
        return new String(Files.readAllBytes(file.toPath()), StandardCharsets.UTF_8);
    }

    private static String tail(final String value) {
        final int start = Math.max(0, value.length() - 4000);
        return value.substring(start);
    }

    private static final class HmcScreenLogRow {
        final long state;
        final double joint;
        final double prior;
        final double likelihood;
        final double rootAge;
        final double ucgdMean;

        HmcScreenLogRow(final long state,
                        final double joint,
                        final double prior,
                        final double likelihood,
                        final double rootAge,
                        final double ucgdMean) {
            this.state = state;
            this.joint = joint;
            this.prior = prior;
            this.likelihood = likelihood;
            this.rootAge = rootAge;
            this.ucgdMean = ucgdMean;
        }
    }
}
