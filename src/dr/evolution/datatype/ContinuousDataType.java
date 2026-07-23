/*
 * ContinuousDataType.java
 *
 * Copyright © 2002-2024 the BEAST Development Team
 * http://beast.community/about
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
 *
 */

package dr.evolution.datatype;

/**
 * @author Andrew Rambaut
 *
 * Continuous data type. This is a place holder to allow mixing of continuous with
 * discrete traits. None of the methods will return anything useful.
 */
public class ContinuousDataType extends DataType {

    public static final String DESCRIPTION = "continuous";
    public static final ContinuousDataType INSTANCE = new ContinuousDataType();

    /**
     * Constructor
     */
    public ContinuousDataType(){
        stateCount = 0;
        ambiguousStateCount = 0;
    }

    @Override
    public char[] getValidChars() {
        return null;
    }

    /**
     * @return the description of the data type
     */
    public String getDescription() {
		return DESCRIPTION;
	}

    public int getType(){
        return CONTINUOUS;
    }

}
