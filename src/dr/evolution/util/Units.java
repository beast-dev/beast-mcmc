/*
 * Units.java
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

package dr.evolution.util;

/**
 * interface holding unit constants
 *
 * @version $Id: Units.java,v 1.17 2005/05/24 20:25:57 rambaut Exp $
 *
 * @author Alexei Drummond
 * @author Andrew Rambaut
 */
public interface Units {
    
	final int SUBSTITUTIONS = 0;
	final int GENERATIONS = 1;
	final int DAYS = 2;
	final int MONTHS = 3;
	final int YEARS = 4;

    // array of unit names
    // second dimension is to allow synonyms -- first element is default
	String[][] UNIT_NAMES = {{"substitutions","mutations"}, {"generations"}, {"days"}, {"months"}, {"years"}};

	/**
	 * Gets the units for this object.
	 */
	int getUnits();

	/**
	 * Sets the units for this object.
	 */
	void setUnits(int units);

    public class Utils {

        public static String getDefaultUnitName(int i) {
            return UNIT_NAMES[i][0];
        }
    }
}

/* might be better to go to something like this:

public interface Units
{
	public static final UnitType EXPECTED_SUBSTITUTIONS = new UnitType("Expected Substitutions");
	public static final UnitType GENERATIONS = new UnitType("Generations");
	public static final UnitType DAYS = new UnitType("Days");
	public static final UnitType MONTHS = new UnitType("Months");
	public static final UnitType YEARS = new UnitType("Years");

	UnitType getUnits();

	void setUnits(UnitType units);

	public class UnitType {
	
		private final String name;
		
		private UnitType(String name) { this.name = name; }
		
		public String toString() { return name; }
	}
}

*/