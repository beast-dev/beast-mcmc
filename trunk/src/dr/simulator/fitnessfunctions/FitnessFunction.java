/*
 * FitnessFunction.java
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
 * An abstract class defining a general fitness function on a byte-encoded molecular sequence
 *
 * @author Alexei Drummond
 * @author Andrew Rambaut
 */
public interface FitnessFunction {

    /**
     * @return the log fitness of the given sequence
     */
	double getLogFitness(byte[] sequence);

	/**
	 * @return the log fitness of the given sequence
	 */
	double getLogFitness(int position, byte state);

    /**
     * @return a sequence which represents the fittest possible
     */
    byte[] getFittestSequence();


}