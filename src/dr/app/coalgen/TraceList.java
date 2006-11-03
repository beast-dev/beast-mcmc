/*
 * TraceList.java
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

package dr.app.coalgen;

/**
 * An interface and default class that stores a set of traces from a single chain
 *
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @version $Id: TraceList.java,v 1.2 2004/10/01 22:40:02 alexei Exp $
 */
 
public interface TraceList {

	/** @return the name of this trace list */
	String getName();

	/** @return the number of traces in this trace list */
	int getTraceCount();
	
	/** @return the index of the trace with the given name */
	int getTraceIndex(String name);
	     
	/** @return the name of the trace with the given index */
	String getTraceName(int index);
	
	/** @return the burn-in for this trace list */
	int getBurnIn();

	/** @return the number of states in the traces */
	int getStateCount();
	
	/** @return the size of the step between states */
	int getStepSize();
	
	/** get the values of trace with the given index (without burnin) */
	void getValues(int index, double[] destination);
}