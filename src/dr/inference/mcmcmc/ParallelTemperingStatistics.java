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

import java.util.ArrayList;

/**
 * @author Marc A. Suchard
 * @author Guy Baele
 */

public class ParallelTemperingStatistics {

    private final StatisticsList[][] allStatistics;
    private final double[] temperatures;

    public ParallelTemperingStatistics(MCMCMCOptions options) {
        double[] temp = options.getChainTemperatures();
        final int length = temp.length;
        this.temperatures = new double[length];
        System.arraycopy(temp, 0, this.temperatures, 0, length); // deep copy to debug

        this.allStatistics = new StatisticsList[length][length];
        for (int i = 0; i < length; ++i) {
            for (int j = 0; j < length; ++j) {
                this.allStatistics[i][j] = new StatisticsList();
            }
        }
    }

    void recordStatistics(int index1, int index2,
                          int rank1, int rank2,
                          double temp1, double temp2,
                          double criterion, boolean success) {
        if (temp1 != temperatures[rank1] || temp2 != temperatures[rank2]) {
            throw new IllegalArgumentException("Temperature mismatch");
        }

        allStatistics[rank1][rank2].add(new Statistics(criterion, success));

        if (DEBUG) {
            System.out.println(rank1 + " (" + index1 + ") <-> " + rank2 + " (" + index2 + ") : " + success);
        }
    }

    public String getReport() {
        return "Parallel tempering report"; // TODO
    }

    static class Statistics {

        double criterion;
        boolean success;

        Statistics(double criterion, boolean success) {
            this.criterion = criterion;
            this.success = success;
        }
    }

    private static final boolean DEBUG = false;

    static class StatisticsList extends ArrayList<Statistics> {}
}
