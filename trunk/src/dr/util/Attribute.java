/*
 * Attribute.java
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

package dr.util;


/**
 * An immutable attribute has a name and value.
 *
 * @author Alexei Drummond
 * @version $Id: Attribute.java,v 1.24 2005/05/24 20:26:01 rambaut Exp $
 */


public interface Attribute {

    public final static String ATTRIBUTE = "att";
    public final static String NAME = "name";
    public final static String VALUE = "value";

    String getAttributeName();

    Object getAttributeValue();

    public class Default implements Attribute {


        public Default(String name, Object value) {
            this.name = name;
            this.value = value;
        }

        public String getAttributeName() {
            return name;
        }

        public Object getAttributeValue() {
            return value;
        }

        public String toString() {
            return name + ": " + value;
        }

        private String name;
        private Object value;

    }

    public class DoubleAttribute extends Default {

        public DoubleAttribute(String name, double value) {
            super(name, java.lang.Double.toString(value));
        }
    }
}

