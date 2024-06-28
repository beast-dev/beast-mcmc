/*
 * TraceTypeUtils.java
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

package dr.inference.trace;

import java.util.List;

/**
 * @author Walter Xie
 */
public class TraceTypeUtils {

    /**
     * return true, if all selected traces are discrete
     * {@link dr.inference.trace.TraceType#isDiscrete() isDiscrete}.
     *
     * @param traceLists
     * @param traceNames
     * @return
     */
    public static boolean allDiscrete(TraceList[] traceLists, List<String> traceNames) {
        for (TraceList traceList : traceLists) {
            for (String traceName : traceNames) {
                int index = traceList.getTraceIndex(traceName);
                TraceDistribution td = traceList.getCorrelationStatistics(index);
                if (td == null) return false;
                if (! td.getTraceType().isDiscrete()) return false;
            }
        }
        return true;
    }

    /**
     * return true, if all selected traces are categorical
     * {@link dr.inference.trace.TraceType#isCategorical() isCategorical}.
     *
     * @param traceLists
     * @param traceNames
     * @return
     */
    public static boolean allCategorical(TraceList[] traceLists, List<String> traceNames) {
        for (TraceList traceList : traceLists) {
            for (String traceName : traceNames) {
                int index = traceList.getTraceIndex(traceName);
                TraceDistribution td = traceList.getCorrelationStatistics(index);
                if (td == null) return false;
                if (! td.getTraceType().isCategorical()) return false;
            }
        }
        return true;
    }

    public static boolean anyCategorical(List<TraceList> traceLists, List<String> traceNames) {
        for (TraceList traceList : traceLists) {
            if (traceNames == null) {
                for (int i = 0; i < traceList.getTraceCount(); i++) {
                    TraceDistribution td = traceList.getCorrelationStatistics(i);
                    if (td != null && ! td.getTraceType().isCategorical())
                        return true;
                }
            } else {
                for (String traceName : traceNames) {
                    int index = traceList.getTraceIndex(traceName);
                    TraceDistribution td = traceList.getCorrelationStatistics(index);
                    if (td != null && ! td.getTraceType().isCategorical())
                        return true;
                }
            }
        }
        return false;
    }


    /**
     * return true, if all selected traces are numbers
     * {@link dr.inference.trace.TraceType#isNumber() isNumber}.
     *
     * @param traceLists
     * @param traceNames
     * @return
     */
    public static boolean allNumeric(TraceList[] traceLists, List<String> traceNames) {
        for (TraceList traceList : traceLists) {
            for (String traceName : traceNames) {
                int index = traceList.getTraceIndex(traceName);
                TraceDistribution td = traceList.getCorrelationStatistics(index);
                if (td == null) return false;
                if (! td.getTraceType().isNumber()) return false;
            }
        }
        return true;
    }


    public static boolean anyNumeric(List<TraceList> traceLists, List<String> traceNames) {
        for (TraceList traceList : traceLists) {
            if (traceNames == null) {
                for (int i = 0; i < traceList.getTraceCount(); i++) {
                    TraceDistribution td = traceList.getCorrelationStatistics(i);
                    if (td != null && ! td.getTraceType().isNumber())
                        return true;
                }
            } else {
                for (String traceName : traceNames) {
                    int index = traceList.getTraceIndex(traceName);
                    TraceDistribution td = traceList.getCorrelationStatistics(index);
                    if (td != null && ! td.getTraceType().isNumber())
                        return true;
                }
            }
        }
        return false;
    }

}
