/*
 * MCMCOptions.java
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

package dr.inference.mcmc;

/**
 * A class that brings together the auxillary information associated
 * with an MCMC analysis.
 *
 * @author Alexei Drummond
 *
 * @version $Id: MCMCOptions.java,v 1.7 2005/05/24 20:25:59 rambaut Exp $
 */
public class MCMCOptions {

	int chainLength;
	int logEvery, repeats = 1;
	boolean verbose;
	boolean coercion = true;
	boolean append = false;
	int preBurnin = 0;
	boolean preBurninSet = false;
	double temperature = 1.0;

	public MCMCOptions() { verbose = true; }
	/** @return the chain length of the MCMC analysis */
	public final int getChainLength() { return chainLength; }

	/** @return whether this MCMC run is verbose. */
	public final boolean isVerbose() {return verbose; }

	public final boolean useCoercion() {return coercion; }

	public final int getRepeats() { return repeats; }
	public final boolean getAppend() { return append; }

	public final int getPreBurnin() { return preBurnin; }

    public final double getTemperature() { return temperature; }
    public final void setTemperature(double temperature) { this.temperature = temperature; }

	public final void setChainLength(int length) {
		chainLength = length;
		if (!preBurninSet) preBurnin=chainLength/100;
	}
	public final void setAppend(boolean append) { this.append = append; }

	public final void setVerbose(boolean ver) { verbose = ver; }
	public final void setRepeats(int reps) { repeats = reps; }

	public final void setUseCoercion(boolean coercion) {
		this.coercion = coercion;
		if (!coercion) preBurnin = 0;
	}

	public final void setPreBurnin(int preBurnin) {
		this.preBurnin = preBurnin;
		preBurninSet = true;
	}
}
