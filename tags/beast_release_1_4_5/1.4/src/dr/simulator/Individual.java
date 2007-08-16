/*
 * Individual.java
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
 *         Date: Apr 22, 2005
 *         Time: 2:23:14 PM
 */
public class Individual {

    public Individual() {
    }

    public Individual(Genome genome, Individual parent) {
        this.genome = genome;
        this.parent = parent;
    }

    public Genome getGenome() {
        return genome;
    }

    public Individual getParent() {
        return parent;
    }

    public void setGenome(Genome genome) {
        this.genome = genome;
    }

    public void setParent(Individual parent) {
        this.parent = parent;
    }

    public int getOffspringCount() {
        return offspringCount;
    }

    public void setOffspringCount(int offspringCount) {
        this.offspringCount = offspringCount;
    }

	public void incrementOffspringCount() {
	    this.offspringCount++;
	}

    public double getLogFitness() {
        return genome.getLogFitness();
    }

    private Genome genome = null;
    private Individual parent = null;
    private int offspringCount = 0;
}
