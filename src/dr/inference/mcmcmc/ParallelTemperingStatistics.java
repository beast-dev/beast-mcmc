/*
 * ParallelTemperingStatistics.java
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

package dr.inference.mcmcmc;

/**
 * Counts attempted and accepted temperature swaps between each pair of ranks.
 *
 * @author Marc A. Suchard
 * @author Guy Baele
 */

public class ParallelTemperingStatistics {

    private final long[][] attempts;
    private final long[][] accepts;
    private final double[][] sumAcceptanceProbability;
    private final double[] temperatures;

    public ParallelTemperingStatistics(MCMCMCOptions options) {
        double[] temp = options.getChainTemperatures();
        final int length = temp.length;
        this.temperatures = new double[length];
        System.arraycopy(temp, 0, this.temperatures, 0, length); // deep copy to debug

        this.attempts = new long[length][length];
        this.accepts = new long[length][length];
        this.sumAcceptanceProbability = new double[length][length];
    }

    void recordStatistics(int index1, int index2,
                          int rank1, int rank2,
                          double temp1, double temp2,
                          double criterion, boolean success) {
        if (temp1 != temperatures[rank1] || temp2 != temperatures[rank2]) {
            throw new IllegalArgumentException("Temperature mismatch");
        }

        int low = Math.min(rank1, rank2);
        int high = Math.max(rank1, rank2);
        attempts[low][high] += 1;
        if (success) {
            accepts[low][high] += 1;
        }
        sumAcceptanceProbability[low][high] += Math.min(1.0, Math.exp(criterion));

        if (DEBUG) {
            System.out.println(rank1 + " (" + index1 + ") <-> " + rank2 + " (" + index2 + ") : " + success);
        }
    }

    /**
     * @return the proportion of accepted swaps over all rank pairs
     */
    public double getAcceptanceRate() {
        long totalAttempts = 0;
        long totalAccepts = 0;
        for (int i = 0; i < attempts.length; ++i) {
            for (int j = 0; j < attempts.length; ++j) {
                totalAttempts += attempts[i][j];
                totalAccepts += accepts[i][j];
            }
        }
        return totalAttempts > 0 ? (double) totalAccepts / totalAttempts : 0.0;
    }

    public String getReport() {
        StringBuilder sb = new StringBuilder();
        sb.append("Parallel tempering swaps:\n");
        sb.append(String.format("  %-10s %-20s %10s %10s %8s %10s%n",
                "ranks", "temperatures", "attempts", "accepted", "rate", "mean P(a)"));
        for (int i = 0; i < attempts.length; ++i) {
            for (int j = i + 1; j < attempts.length; ++j) {
                long n = attempts[i][j];
                if (n > 0) {
                    sb.append(String.format("  %-10s %-20s %10d %10d %8.4f %10.4f%n",
                            i + " <-> " + j,
                            String.format("%.4f <-> %.4f", temperatures[i], temperatures[j]),
                            n, accepts[i][j], (double) accepts[i][j] / n,
                            sumAcceptanceProbability[i][j] / n));
                }
            }
        }
        return sb.toString();
    }

    private static final boolean DEBUG = false;
}
