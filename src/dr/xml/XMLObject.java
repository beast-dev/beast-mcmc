/*
 * XMLObject.java
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

package dr.xml;

import dr.inference.model.AbstractModel;
import org.w3c.dom.Element;
import org.w3c.dom.NamedNodeMap;

import java.lang.reflect.Constructor;
import java.util.ArrayList;
import java.util.List;
import java.util.StringTokenizer;
import java.util.Vector;

/**
 * This class wraps a DOM Element for the purposes of parsing.
 *
 * @author Alexei Drummond
 * @version $Id: XMLObject.java,v 1.30 2005/07/11 14:06:25 rambaut Exp $
 */
public class XMLObject {

    public static final String missingValue = "NA";

    /**
     * @param e the element the construct this XML object from
     */
    public XMLObject(Element e, XMLObject parent) {
        this.element = e;
        this.parent = parent;
    }

    public XMLObject(XMLObject obj, int index) {

       this(obj.element, null);
       nativeObject = ((List)obj.getNativeObject()).get(index);
   }

    /**
     * @return the number of children this XMLObject has.
     */
    public final int getChildCount() {
        return children.size();
    }

    /**
     * @param i the index of the child to return
     * @return the ith child in native format if available, otherwise as
     *         an XMLObject.
     */
    public Object getChild(int i) {
        Object obj = getRawChild(i);
        XMLObject xo = null;

        if( obj instanceof XMLObject ) {
            xo = (XMLObject) obj;
        } else if( obj instanceof Reference ) {
            xo = ((Reference) obj).getReferenceObject();
        }

        if( xo != null && xo.hasNativeObject() ) {
            return xo.getNativeObject();
        }
        return obj;
    }

    /**
     * @param c the class of the child to return
     * @return the first child with a native format of the given class, or null if no such child exists.
     */
    public Object getChild(Class c) {

        for (int i = 0; i < getChildCount(); i++) {
            Object child = getChild(i);
            if( c.isInstance(child) ) {
                return child;
            }
        }
        return null;
    }

    /**
     * @return all children with or empty list if no children.
     */
    public List<Object> getChildren() {

        List<Object> allChildren = new ArrayList<Object>();
        for (int i = 0; i < getChildCount(); i++) {
            Object child = getChild(i);
                allChildren.add(child);

        }
        return allChildren;
    }

    /**
     * @param c the class of the children to return
     * @return all children with a native format of the given class, or empty if no such child exists.
     */
    public <T> List<T> getAllChildren(Class<T> c) {

        List<T> allChildren = new ArrayList<T>();;
        for (int i = 0; i < getChildCount(); i++) {
            Object child = getChild(i);
            if( c.isInstance(child) ) {
                allChildren.add(c.cast(child));
            }

        }
        return allChildren;

    }

    /**
     * @param name the name of the child to return
     * @return the first child of type XMLObject with a given name, or null if no such child exists.
     */
    public XMLObject getChild(String name) {

        for (int i = 0; i < getChildCount(); i++) {
            Object child = getChild(i);
            if( child instanceof XMLObject ) {
                if( ((XMLObject) child).getName().equals(name) ) {
                    return (XMLObject) child;
                }
            }
        }
        return null;
    }

    /**
     * @param name the name of the children
     * @return all children with a given name.
     */
    public List<XMLObject> getAllChildren(String name) {

        List<XMLObject> allChildren = new ArrayList<XMLObject>();
        for (int i = 0; i < getChildCount(); i++) {
            Object child = getChild(i);
            if( child instanceof XMLObject && ((XMLObject)child).getName().equals(name)) {
                allChildren.add((XMLObject)child);
            }
        }
        return allChildren;

    }

    /**
     * @param elementName the name of the XML wrapper element the child resides in
     * @return the first child element out of a named XMLObject element
     * @throws XMLParseException if no wrapper element exists, or it the child in
     *                           wrapper element is not an XMLObject
     */
    public Object getElementFirstChild(String elementName) throws XMLParseException {
        Object child = getChild(elementName);

        if (child == null)
            throw new XMLParseException("Child element called " + elementName +
                    " does not exist inside element " + getName());
        if (!(child instanceof XMLObject))
            throw new XMLParseException("Child element called " + elementName +
                    " inside element " + getName() + " is not an XMLObject.");

        return ((XMLObject) child).getChild(0);
    }

    public String getChildName(int i) {

        Object obj = getRawChild(i);
        XMLObject xo;

        if (obj instanceof XMLObject) {
            xo = (XMLObject) obj;
        } else if (obj instanceof Reference) {
            xo = ((Reference) obj).getReferenceObject();
        } else {
            return "";
        }

        return xo.getName();
    }

    /**
     * @param name the name of the child element being tested for.
     * @return true if a child element of the given name exists.
     */
    public boolean hasChildNamed(String name) {
        final Object child = getChild(name);
        return (child != null) && (child instanceof XMLObject);
    }

    /**
     * @return all attributes
     */
    public NamedNodeMap getAttributes() {
        return element.getAttributes();
    }

    /**
     * @param i the index of the child to return
     * @return the ith child as a boolean.
     * @throws XMLParseException if getChild(i) would
     */
    public boolean getBooleanChild(int i) throws XMLParseException {
        return getBoolean(getChild(i));
    }

    /**
     * @param i the index of the child to return
     * @return the ith child as a double.
     * @throws XMLParseException if getChild(i) would
     */
    public double getDoubleChild(int i) throws XMLParseException {
        return getDouble(getChild(i));
    }

    /**
     * @param i the index of the child to return
     * @return the ith child as a double[].
     * @throws XMLParseException if getChild(i) would
     */
    public double[] getDoubleArrayChild(int i) throws XMLParseException {
        return getDoubleArray(getChild(i));
    }

    /**
     * @param i the index of the child to return
     * @return the ith child as an integer.
     * @throws XMLParseException if getChild(i) would
     */
    public int getIntegerChild(int i) throws XMLParseException {
        return getInteger(getChild(i));
    }

    /**
     * @param i the index of the child to return
     * @return the ith child as a string.
     * @throws XMLParseException if getChild(i) would
     */
    public String getStringChild(int i) throws XMLParseException {
        return getString(getChild(i));
    }

    /**
     * @param i the index of the child to return
     * @return the ith child as a String[].
     * @throws XMLParseException if getChild(i) would
     */
    public String[] getStringArrayChild(int i) throws XMLParseException {
        return getStringArray(getChild(i));
    }

    /**
     * Attribute value, if present - default otherwise.
     *
     * @param name         attribute name
     * @param defaultValue the default value
     * @return the given attribute if it exists, otherwise the default value
     * @throws XMLParseException if attribute can't be converted to desired type
     */
    public <T> T getAttribute(String name, T defaultValue) throws XMLParseException {
        if (element.hasAttribute(name)) {

            final String s = element.getAttribute(name);
            for (Constructor c : defaultValue.getClass().getConstructors()) {
                final Class[] classes = c.getParameterTypes();
                if (classes.length == 1 && classes[0].equals(String.class)) {
                    try {
                        return (T) c.newInstance(s);
                    } catch (Exception e) {
                        throw new XMLParseException(" conversion of '" + s + "' to " +
                                defaultValue.getClass().getName() + " failed");
                    }
                }
            }
        }
        return defaultValue;
    }

    /**
     * @return the named attribute
     */
    public Object getAttribute(String name) throws XMLParseException {
        return getAndTest(name);
    }

    /**
     * @return the named attribute as a boolean.
     */
    public boolean getBooleanAttribute(String name) throws XMLParseException {
        return getBoolean(getAndTest(name));
    }

    /**
     * @return the named attribute as a double.
     */
    public double getDoubleAttribute(String name) throws XMLParseException {
        return getDouble(getAndTest(name));
    }

    /**
     * @return the named attribute as a double[].
     */
    public double[] getDoubleArrayAttribute(String name) throws XMLParseException {
        return getDoubleArray(getAndTest(name));
    }

    /**
     * @return the named attribute as a double[].
     */
    public int[] getIntegerArrayAttribute(String name) throws XMLParseException {
        return getIntegerArray(getAndTest(name));
    }

    /**
     * @return the named attribute as an integer.
     */
    public int getIntegerAttribute(String name) throws XMLParseException {
        return getInteger(getAndTest(name));
    }

    /**
     * @return the named attribute as a long integer.
     */
    public long getLongIntegerAttribute(String name) throws XMLParseException {
        return getLongInteger(getAndTest(name));
    }

    /**
     * @return the named attribute as a string.
     */
    public String getStringAttribute(String name) throws XMLParseException {
        return getString(getAndTest(name));
    }

    /**
     * @return the named attribute as a String[].
     */
    public String[] getStringArrayAttribute(String name) throws XMLParseException {
        return getStringArray(getAndTest(name));
    }

    /**
     * @param valueList if this ArrayList is not null it is populated by the double array
     *                  that the given string represents.
     * @return true if the given string represents a whitespaced-delimited array of doubles.
     */
    public static boolean isDoubleArray(String s, List<Double> valueList) {
        try {
            StringTokenizer st = new StringTokenizer(s);
            while (st.hasMoreTokens()) {
                String token = st.nextToken();
                Double d;
                if (token.compareToIgnoreCase(missingValue) == 0)
                    d = Double.NaN;
                else
                    d = new Double(token);
                if (valueList != null) valueList.add(d);
            }
            return true;
        } catch (NumberFormatException e) {
            return false;
        }
    }

    /**
     * @param valueList if this ArrayList is not null it is populated by the integer array
     *                  that the given string represents.
     * @return true if the given string represents a whitespaced-delimited array of integers.
     */
    public static boolean isIntegerArray(String s, List<Integer> valueList) {
        try {
            StringTokenizer st = new StringTokenizer(s);
            while (st.hasMoreTokens()) {
                Integer d = new Integer(st.nextToken());
                if (valueList != null) valueList.add(d);
            }
            return true;
        } catch (NumberFormatException e) {
            return false;
        }
    }

    public final static String ID = "id";

    public boolean hasId() {
        return hasAttribute(ID);
    }

    public String getId() throws XMLParseException {
        return getStringAttribute(ID);
    }

    /**
     * @return true if either an attribute exists.
     */
    public boolean hasAttribute(String name) {
        return (element.hasAttribute(name));
    }

    public String getName() {
        return element.getTagName();
    }

    public Object getNativeObject() {
        return nativeObject;
    }

    public boolean hasNativeObject() {
        return nativeObject != null;
    }

    public String toString() {

        String prefix = getName();
        if (hasId()) {
            try {
                prefix += ":" + getId();
            } catch (XMLParseException e) {
                // this shouldn't happen
                assert false;
            }
        }
        //if (nativeObject != null) return prefix + ":" + nativeObject.toString();
        return prefix;
    }

    public String content() {

        if (nativeObject != null) {
            if (nativeObject instanceof dr.util.XHTMLable) {
                return ((dr.util.XHTMLable) nativeObject).toXHTML();
            } else {
                return nativeObject.toString();
            }
        }
        return "";
    }

    //*********************************************************************
    // Package functions
    //*********************************************************************

    /**
     * Adds a child.
     */
    void addChild(Object child) {
        if (child instanceof XMLObject ||
                child instanceof Reference ||
                child instanceof String) {

            children.add(child);
        } else throw new IllegalArgumentException();
    }

    /**
     * @return the ith child of this XMLObject, without processing.
     */
    public Object getRawChild(int i) {
        return children.get(i);
    }

    /**
     * Sets the native object represented by this XMLObject.
     */
    public void setNativeObject(Object obj) {
        nativeObject = obj;
    }

    boolean isReference(int child) {
        return (getRawChild(child) instanceof Reference);
    }

    //*********************************************************************
    // Static members
    //*********************************************************************

    //*********************************************************************
    // Private methods
    //*********************************************************************

    /**
     * @return the object as a boolean if possible
     */
    private boolean getBoolean(Object obj) throws XMLParseException {

        if (obj instanceof Boolean) return (Boolean) obj;
        if (obj instanceof String) {
            if (obj.equals("true")) return true;
            if (obj.equals("false")) return false;
        }
        throw new XMLParseException("Expected a boolean (true|false), but got " + obj);
    }

    /**
     * @return the object as an double if possible
     */
    private double getDouble(Object obj) throws XMLParseException {
        try {
            if (obj instanceof Number) {
                return ((Number) obj).doubleValue();
            }
            if (obj instanceof String) {
                return Double.parseDouble((String) obj);
            }
        } catch (NumberFormatException nfe) {
            throw new XMLParseException("Expected double precision number, but got " + obj);
        }
        throw new XMLParseException("Expected double precision number, but got " + obj);
    }

    /**
     * @return the object as an double[] if possible
     */
    private double[] getDoubleArray(Object obj) throws XMLParseException {

        if (obj instanceof Number) return new double[]{((Number) obj).doubleValue()};
        if (obj instanceof double[]) return (double[]) obj;
        if (obj instanceof String) {
            List<Double> valueList = new ArrayList<Double>();
            if (isDoubleArray((String) obj, valueList)) {
                double[] values = new double[valueList.size()];
                for (int i = 0; i < values.length; i++) {
                    values[i] = valueList.get(i);
                }
                return values;
            } else {
                throw new XMLParseException("Expected array of double precision numbers, but got " + obj);
            }
        }
        throw new XMLParseException("Expected array of double precision numbers, but got " + obj);
    }

    /**
     * @return the object as an double[] if possible
     */
    private int[] getIntegerArray(Object obj) throws XMLParseException {

        if (obj instanceof Number) return new int[]{((Number) obj).intValue()};
        if (obj instanceof int[]) return (int[]) obj;
        if (obj instanceof String) {
            ArrayList<Integer> valueList = new ArrayList<Integer>();
            if (isIntegerArray((String) obj, valueList)) {
                int[] values = new int[valueList.size()];
                for (int i = 0; i < values.length; i++) {
                    values[i] = valueList.get(i);
                }
                return values;
            } else {
                throw new XMLParseException("Expected array of integers, but got " + obj);
            }
        }
        throw new XMLParseException("Expected array of integers, but got " + obj);
    }

    /**
     * @return the object as an integer if possible
     */
    private int getInteger(Object obj) throws XMLParseException {
        if (obj instanceof Number) return ((Number) obj).intValue();
        try {
            return Integer.parseInt((String) obj);
        } catch (NumberFormatException e) {
            throw new XMLParseException("Expected integer, got " + obj);
        }
    }

    /**
     * @return the object as an integer if possible
     */
    private long getLongInteger(Object obj) throws XMLParseException {
        if (obj instanceof Number) return ((Number) obj).longValue();
        try {
            return Long.parseLong((String) obj);
        } catch (NumberFormatException e) {
            throw new XMLParseException("Expected long integer, got " + obj);
        }
    }

    /**
     * @return the object as a string if possible
     */
    private String getString(Object obj) throws XMLParseException {

        if (obj instanceof String) return (String) obj;
        throw new XMLParseException("Expected string, but got " + obj);
    }

    /**
     * @return the object as an double[] if possible
     */
    private String[] getStringArray(Object obj) throws XMLParseException {

        if (obj instanceof String[]) return (String[]) obj;
        if (obj instanceof String) {
            ArrayList<String> stringList = new ArrayList<String>();
            StringTokenizer st = new StringTokenizer((String) obj);
            while (st.hasMoreTokens()) {
                stringList.add(st.nextToken());
            }
            String[] strings = new String[stringList.size()];
            for (int i = 0; i < strings.length; i++) {
                strings[i] = stringList.get(i);
            }
            return strings;
        }
        throw new XMLParseException("Expected array of strings, but got " + obj);
    }

    /**
     * @return the named attribute if it exists, throws XMLParseException otherwise.
     */
    private Object getAndTest(String name) throws XMLParseException {

        if (element.hasAttribute(name)) {
            return element.getAttribute(name);
        }
        throw new XMLParseException("'" + name + "' attribute was not found in " + element.getTagName() + " element.");
    }

    public XMLObject getParent() {
        return parent;
    }

    //*********************************************************************
    // Private instance variables
    //*********************************************************************

    private final Vector<Object> children = new Vector<Object>();
    private final Element element;
    private final XMLObject parent;

    private Object nativeObject;

    // The objectStore representing the local scope of this element.
//	private ObjectStore store;
}
