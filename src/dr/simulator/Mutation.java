/*
 * Mutation.java
 *
 * Copyright (C) 2002-2006 Alexei Drummond and Andrew Rambaut
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

package dr.simulator;

/**
 * @author rambaut
 *         Date: Apr 26, 2005
 *         Time: 10:28:35 AM
 */
public final class Mutation {
    private Mutation(int position, byte state) {
        this.position = position;
        this.state = state;
    }

    final int position;
    final byte state;

    public static Mutation getMutation(int position, byte state) {
        return mutations[position][state];
    }

    private static Mutation[][] mutations = null;

    public static void initialize(int genomeLength, int stateSize) {
        mutations = new Mutation[genomeLength][];
        for (int i = 0; i < genomeLength; i++) {
            mutations[i] = new Mutation[stateSize];
            for (byte j = 0; j < stateSize; j++) {
                mutations[i][j] = new Mutation(i, j);
            }
        }
    }

}
