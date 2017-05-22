/*
 * TraceFactory.java
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

package dr.inference.trace;

/**
 * The trace type to determine what statistics
 *
 * @author Alexei Drummond
 * @author Walter Xie
 * @author Andrew Rambaut
 */
public enum TraceType {
    // changed this to 'real' as this is less Comp Sci. than 'double'
    REAL("real", "R", Double.class),
    INTEGER("integer", "I", Integer.class),
    CATEGORICAL("categorical", "C", String.class),
    BINARY("binary", "B", Integer.class);

    TraceType(String name, String brief, Class type) {
        this.name = name;
        this.brief = brief;
        this.type = type;
    }

    public String toString() {
        return name;
    }

    public String getBrief() {
        return brief;
    }

    private Class getTypeClass() {
        return type;
    }

    public boolean isContinuous() {
        return getTypeClass() == Double.class;
    }

    /**
     * is the data type numerical
     * @return
     */
    public boolean isNumber() {
        return getTypeClass() == Double.class || getTypeClass() == Integer.class;
    }

    /**
     * is the datatype ordinal (integer or binary)
     * @return
     */
    public boolean isInteger() {
        return getTypeClass() == Integer.class;
    }

    /**
     * is the datatype discrete (ordinal or categorical)
     * @return
     */
    public boolean isDiscrete() {
        return isInteger() || isCategorical();
    }

    /**
     * is the datatype binary
     * @return
     */
    public boolean isBinary() {
        return this == BINARY;
    }

    public boolean isIntegerOrBinary() {
        return isInteger() || isBinary();
    }

    /**
     * is the datatype categorical
     * @return
     */
    public boolean isCategorical() {
        return getTypeClass() == String.class;
    }

    private final String name;
    private final String brief;
    private final Class type;

}
