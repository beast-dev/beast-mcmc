/*
 * HistoryFilter.java
 *
 * Copyright (c) 2002-2015 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.evomodel.treelikelihood.utilities;

import java.util.Set;

/**
 * @author Marc A. Suchard
 */
public interface HistoryFilter {
    public boolean filter(String source, String destination, double time);

    public String getDescription();

    public class Default implements HistoryFilter {
        public Default() {
            // Do nothing
        }

        public String getDescription() {
            return "Default";
        }

        public boolean filter(String source, String destination, double time) {
            return true;
        }
    }

    public class SetFilter implements HistoryFilter {
        final private Set<String> sources;
        final private Set<String> destinations;
        final double maxTime;
        final double minTime;

        public SetFilter(Set<String> sources, Set<String> destinations, double maxTime, double minTime) {
            this.sources = sources;
            this.destinations = destinations;
            this.maxTime = maxTime;
            this.minTime = minTime;
        }

        public String getDescription() {
            StringBuilder sb = new StringBuilder();
            sb.append(minTime).append(" <= event time <= ").append(maxTime);
            return sb.toString();
        }

        public boolean filter(String source, String destination, double time) {

            if (sources == null && destinations == null) {
                return time <= maxTime && time >= minTime;
            } // else

            return sources.contains(source) && destinations.contains(destination) &&
                    time <= maxTime && time >= minTime;
        }
    }
}
