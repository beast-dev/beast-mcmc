/*
 * Variable.java
 *
 * Copyright (C) 2002-2009 Alexei Drummond and Andrew Rambaut
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
 * BEAST is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with BEAST; if not, write to the
 * Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
 * Boston, MA  02110-1301  USA
 */

package dr.inference.model;

/**
 * A generic random variable.
 *
 * @author Alexei Drummond
 */
public interface Variable<V> {

    public enum ChangeType {
        VALUE_CHANGED,
        REMOVED,
        ADDED
    }

    /**
     * @return the name of this variable.
     */
    public String getVariableName();

    public V getValue(int index);

    public void setValue(int index, V value);

    /**
     * @return the size of this variable - i.e. the length of the vector
     */
    public int getSize();

    /**
     * adds a parameter listener that is notified when this parameter changes.
     *
     * @param listener the listener
     */
    void addVariableListener(VariableListener listener);

    /**
     * removes a parameter listener.
     *
     * @param listener the listener
     */
    void removeVariableListener(VariableListener listener);

    /**
     * stores the state of this parameter for subsquent restore
     */
    void storeVariableValues();

    /**
     * restores the stored state of this parameter
     */
    void restoreVariableValues();

    /**
     * accepts the stored state of this parameter
     */
    void acceptVariableValues();
}
