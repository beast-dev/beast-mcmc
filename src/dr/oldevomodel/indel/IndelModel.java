/*
 * IndelModel.java
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

package dr.oldevomodel.indel;

import dr.evolution.util.Units;
import dr.inference.model.AbstractModel;



/**
 * This abstract class contains methods that are of general use for
 * modelling insertions and deletions.
 *
 * @version $$
 *
 * @author Alexei Drummond
 */
public abstract class IndelModel extends AbstractModel implements Units {

	public static final String INDEL_MODEL = "indelModel";

	//
	// Public stuff
	//

	public IndelModel(String name) { 
		super(name);
		
		units = Units.Type.GENERATIONS;
	}

	//
	// functions that define an indel model (left for subclass)
	//

	/**
	 * Gets the birth rate of insertions of a given length.
	 */
	public abstract double getBirthRate(int length);

	/**
	 * Gets the death rate of insertions of a given length.
	 */
	public abstract double getDeathRate(int length);
	
	/**
	 * Units in which time units are measured.
	 */
	private Type units;

	/**
	 * sets units of measurement.
	 *
	 * @param u units
	 */
	public void setUnits(Type u)
	{
		units = u;
	}

	/**
	 * returns units of measurement.
	 */
	public Type getUnits()
	{
		return units;
	}
}
