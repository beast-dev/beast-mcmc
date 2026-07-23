/*
 * IntervalChangedEvent.java
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

package dr.evomodel.bigfasttree;

/**
 * @author Marc A. Suchard
 * @author Guy Baele
 */
public interface IntervalChangedEvent {

    int getFirstInterval();
    int getLastInterval();

    class FirstAffectedInterval implements IntervalChangedEvent {

        private final int firstInterval;

        public FirstAffectedInterval(int interval) {
            this.firstInterval = interval;
        }

        @Override
        public int getFirstInterval() {
            return firstInterval;
        }

        @Override
        public int getLastInterval() {
            throw new RuntimeException("");
        }

    }

    class AffectedIntervals implements IntervalChangedEvent {

        private final int firstInterval;
        private final int lastInterval;

        public AffectedIntervals(int firstInterval, int lastInterval) {
            this.firstInterval = firstInterval;
            this.lastInterval = lastInterval;
        }

        public AffectedIntervals(int[] intervals) {
            this.firstInterval = intervals[0];
            this.lastInterval = intervals[1];
        }

        @Override
        public int getFirstInterval() {
            return firstInterval;
        }

        @Override
        public int getLastInterval() {
            return lastInterval;
        }

    }

}
