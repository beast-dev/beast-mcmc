/*
 * AbstractFitnessFunction.java
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

package dr.simulator.fitnessfunctions;

/**
 * @author rambaut
 *         Date: Apr 28, 2005
 *         Time: 12:24:46 PM
 */
public class AbstractFitnessFunction implements FitnessFunction {

    protected void setLogFitnessTable(double[][] logFitnessTable) {
        this.logFitnessTable = logFitnessTable;

        fittestSequence = new byte[logFitnessTable.length];

        for (int i = 0; i < fittestSequence.length; i++) {
            fittestSequence[i] = 0;
            for (byte j = 1; j < logFitnessTable[0].length; j++) {
                if (logFitnessTable[i][j] > logFitnessTable[i][fittestSequence[i]]) {
                    fittestSequence[i] = j;
                }
            }
        }
    }

    /**
     * @return the log fitness of the given sequence
     */
    public double getLogFitness(byte[] sequence) {
        double totalFitness = 0.0;

        for (int i = 0; i < sequence.length; i++) {
            totalFitness += logFitnessTable[i][sequence[i]];
        }
        return totalFitness;
    }

    /**
     * @return the log fitness of the given sequence
     */
    public double getLogFitness(int position, byte state) {
        return logFitnessTable[position][state];
    }

    /**
     * @return a sequence which represents the fittest possible
     */
    public byte[] getFittestSequence() {
        return fittestSequence;
    }

    private double[][] logFitnessTable;
    private byte[] fittestSequence;
}
