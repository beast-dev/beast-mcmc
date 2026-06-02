/*
 * BenchmarkTimer.java
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

package dr.util;

import java.util.HashMap;
import java.util.Map;

/**
 * @author Marc A. Suchard
 */
public class BenchmarkTimer {

    public void startTimer(String key) {
        startTimes.put(key, getTime());
    }

    public void stopTimer(String key) {

        long end = getTime();
        long start = startTimes.get(key);

        Long total = totals.get(key);
        if (total == null) {
            total = 0L;
        }

        totals.put(key, total + (end - start));
    }
    
    public String toString() {

        StringBuilder sb = new StringBuilder();

        sb.append("\nTIMING:");
        for (String key : totals.keySet()) {
            String value = String.format("%4.3e", (double) totals.get(key));
            sb.append("\n").append(key).append(delimiter).append(value);
        }
        sb.append("\n");

        return sb.toString();
    }

    private long getTime() {
//        return System.nanoTime();
        return System.currentTimeMillis();
    }

    private final Map<String, Long> startTimes = new HashMap<>();
    private final Map<String, Long> totals = new HashMap<>();

    private final static String delimiter = "\t\t";
}
