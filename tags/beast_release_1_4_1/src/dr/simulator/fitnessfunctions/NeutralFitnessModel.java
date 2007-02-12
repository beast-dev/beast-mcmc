/*
 * NeutralFitnessModel.java
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
 * A fitness function in which all molecular sequences have equal and optimal fitness.
 *
 * @author Alexei Drummond
 */
public class NeutralFitnessModel implements FitnessFunction {

    public NeutralFitnessModel(int genomeLength) {
        this.genomeLength = genomeLength;
    }

    /**
     * @return the log fitness of the given sequence
     */
    public double getLogFitness(byte[] sequence) {
        return 0;
    }

    /**
     * @return the log fitness of the given sequence
     */
    public double getLogFitness(int position, byte state) {
        return 0;
    }

    /**
     * @return a sequence which represents the fittest possible
     */
    public byte[] getFittestSequence() {
        // will be initialized to state 0
        return new byte[genomeLength];
    }

    private int genomeLength;
}