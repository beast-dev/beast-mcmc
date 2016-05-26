/*
 * Property.java
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

package dr.util;

import java.lang.reflect.Method;


/**
 * Gets a property of another object using introspection.
 *
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @version $Id: Property.java,v 1.19 2005/05/24 20:26:01 rambaut Exp $
 */
public class Property implements Attribute {

    private Object object = null;
    private Method getter = null;
    private Object argument = null;
    private String name = null;

    public Property(Object object, String name) {
        this(object, name, null);
    }

    public Property(Object object, String name, Object argument) {

        this.name = name;
        this.argument = argument;

        this.object = object;

        StringBuffer getterName = new StringBuffer("get");
        getterName.append(name.substring(0, 1).toUpperCase());
        getterName.append(name.substring(1));
        Class c = object.getClass();

        //System.out.println(getterName + "(" + argument + ")");

        try {
            if (argument != null)
                getter = c.getMethod(getterName.toString(), new Class[]{argument.getClass()});
            else
                getter = c.getMethod(getterName.toString(), (Class[]) null);
        } catch (NoSuchMethodException e) {
        }

    }

    public Method getGetter() {
        return getter;
    }

    //public Object getObject() { return object; }

    public String getAttributeName() {
        if (argument == null) return name;
        return name + "." + argument;
    }

    public Object getAttributeValue() {
        if (object == null || getter == null)
            return null;

        Object result = null;
        Object[] args = null;
        if (argument != null)
            args = new Object[]{argument};

        try {
            result = getter.invoke(object, args);
        } catch (Exception e) {

            e.printStackTrace(System.out);

            throw new RuntimeException(e.getMessage());
        }

        return result;
    }

    public String getPropertyName() {
        return name;
    }

    public String toString() {
        return getAttributeValue().toString();
    }
}
