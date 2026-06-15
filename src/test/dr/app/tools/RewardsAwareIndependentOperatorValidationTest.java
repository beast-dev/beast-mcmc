/*
 * RewardsAwareIndependentOperatorValidationTest.java
 *
 * Copyright (c) 2002-2026 the BEAST Development Team
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
 * BEAST is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with BEAST; if not, write to the
 * Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
 * Boston, MA  02110-1301  USA
 */

package test.dr.app.tools;

import dr.app.tools.RewardsAwareIndependentOperatorValidation;
import dr.inference.operators.RewardsMixtureBranchResamplingHelper;
import test.dr.math.MathTestCase;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;

/**
 * Regression coverage for independent reward-mixture validation scenarios.
 *
 * @author Filippo Monti
 */
public class RewardsAwareIndependentOperatorValidationTest extends MathTestCase {

    public void testContinuousRewardOpenSupportRejectsEndpoints() {
        assertTrue(RewardsMixtureBranchResamplingHelper.isContinuousRewardOutsideOpenSupport(0.0, 0.0, 1.0));
        assertTrue(RewardsMixtureBranchResamplingHelper.isContinuousRewardOutsideOpenSupport(1.0, 0.0, 1.0));
        assertTrue(RewardsMixtureBranchResamplingHelper.isContinuousRewardOutsideOpenSupport(-0.01, 0.0, 1.0));
        assertTrue(RewardsMixtureBranchResamplingHelper.isContinuousRewardOutsideOpenSupport(1.01, 0.0, 1.0));
        assertFalse(RewardsMixtureBranchResamplingHelper.isContinuousRewardOutsideOpenSupport(0.5, 0.0, 1.0));
    }

    public void testScenarioFiveBoundaryProposalsAreRejected() throws Exception {
        final File outputDir = createTempDirectory("reward-validation-s5-boundary");
        try {
            RewardsAwareIndependentOperatorValidation.main(new String[]{
                    outputDir.getAbsolutePath(),
                    "50",
                    "3",
                    "5",
                    "7000",
                    "0",
                    "1000",
                    "777005"
            });

            final double finalLogLikelihood = readSummaryValue(outputDir, "finalLogLikelihood");
            assertTrue("Final likelihood should remain finite after boundary proposals are rejected",
                    Double.isFinite(finalLogLikelihood));
        } finally {
            deleteRecursively(outputDir);
        }
    }

    public void testValidationLogIncludesBurnInRows() throws Exception {
        final File outputDir = createTempDirectory("reward-validation-burnin-log");
        try {
            RewardsAwareIndependentOperatorValidation.main(new String[]{
                    outputDir.getAbsolutePath(),
                    "10",
                    "3",
                    "1",
                    "3000",
                    "2000",
                    "1000",
                    "777101"
            });

            assertEquals("Log should include rows from sweep zero through the final sweep",
                    4, countDataRows(new File(outputDir, "independent_operator_validation.log")));
            assertEquals("Logged sample count should include burn-in rows",
                    4, (int) readSummaryValue(outputDir, "sampleCount"));
        } finally {
            deleteRecursively(outputDir);
        }
    }

    private static File createTempDirectory(final String prefix) throws IOException {
        final File tmp = File.createTempFile(prefix, "");
        if (!tmp.delete()) {
            throw new IOException("Could not delete temporary file: " + tmp);
        }
        if (!tmp.mkdirs()) {
            throw new IOException("Could not create temporary directory: " + tmp);
        }
        return tmp;
    }

    private static double readSummaryValue(final File outputDir,
                                           final String metric) throws IOException {
        final File summary = new File(outputDir, "summary.tsv");
        try (BufferedReader reader = new BufferedReader(new FileReader(summary))) {
            String line;
            while ((line = reader.readLine()) != null) {
                final String[] fields = line.split("\\t");
                if (fields.length == 2 && metric.equals(fields[0])) {
                    return Double.parseDouble(fields[1]);
                }
            }
        }
        throw new IOException("Metric not found in " + summary + ": " + metric);
    }

    private static int countDataRows(final File logFile) throws IOException {
        try (BufferedReader reader = new BufferedReader(new FileReader(logFile))) {
            int rows = 0;
            String line;
            boolean header = true;
            while ((line = reader.readLine()) != null) {
                if (header) {
                    header = false;
                } else if (line.trim().length() > 0) {
                    rows++;
                }
            }
            return rows;
        }
    }

    private static void deleteRecursively(final File file) throws IOException {
        if (file == null || !file.exists()) {
            return;
        }
        if (file.isDirectory()) {
            final File[] children = file.listFiles();
            if (children != null) {
                for (File child : children) {
                    deleteRecursively(child);
                }
            }
        }
        if (!file.delete()) {
            throw new IOException("Could not delete " + file);
        }
    }
}
