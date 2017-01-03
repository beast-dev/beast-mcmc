/*
 * HiddenDataType.java
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

package dr.evolution.datatype;

/**
 * @author Marc Suchard
 */

public interface HiddenDataType {

    int getHiddenClassCount();

    int getStateCount();

    class Utils {

        static boolean[] getStateSet(final int state,
                                     final int stateCount,
                                     final int hiddenClassCount,
                                     final DataType baseDataType) {

//            System.err.println("state = ");
            
            final boolean[] originalStateSet = baseDataType.getStateSet(state);
            boolean[] stateSet = new boolean[stateCount * hiddenClassCount];

            int offset = 0;
            for (int h = 0; h < hiddenClassCount; ++h) {
                System.arraycopy(originalStateSet, 0, stateSet, offset, stateCount);
                offset += stateCount;
            }

            return stateSet;
        }
    }
}
