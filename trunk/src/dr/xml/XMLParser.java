/*
 * XMLParser.java
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

package dr.xml;

import dr.evolution.util.Units;
import dr.util.Identifiable;
import org.w3c.dom.*;
import org.xml.sax.InputSource;

import java.io.Reader;
import java.util.*;

public class XMLParser {

    public static final String ID = "id";
    public static final String IDREF = "idref";
    public static final String CONCURRENT = "concurrent";

    private Vector<Thread> threads = new Vector<Thread>();

    public XMLParser() {
        addXMLObjectParser(new ArrayParser());
        addXMLObjectParser(Report.PARSER);
    }

    public XMLParser(boolean verbose) {
        this();
        this.verbose = verbose;
    }

    public void addXMLObjectParser(XMLObjectParser parser) {

        String[] parserNames = parser.getParserNames();

        for (String parserName : parserNames) {
            XMLObjectParser oldParser = parserStore.get(parserName);
            if (oldParser != null) {
                throw new IllegalArgumentException("New parser ( " + parser + " ) cannot replace existing parser (" + oldParser + ")");
            }

            parserStore.put(parserName, parser);
        }
    }

    public Iterator getParserNames() {
        return parserStore.keySet().iterator();
    }

    public XMLObjectParser getParser(String name) {
        return parserStore.get(name);
    }

    public Iterator getParsers() {
        return parserStore.values().iterator();
    }

    public Iterator getThreads() {
        return threads.iterator();
    }

    public void storeObject(String name, Object object) {

        XMLObject xo = new XMLObject(null /*, objectStore*/);
        xo.setNativeObject(object);
        store.put(name, xo);
    }

    public ObjectStore parse(Reader reader, boolean run)
            throws java.io.IOException,
            org.xml.sax.SAXException,
            dr.xml.XMLParseException,
            javax.xml.parsers.ParserConfigurationException {

        InputSource in = new InputSource(reader);
        javax.xml.parsers.DocumentBuilderFactory documentBuilderFactory = javax.xml.parsers.DocumentBuilderFactory.newInstance();

        javax.xml.parsers.DocumentBuilder documentBuilder = documentBuilderFactory.newDocumentBuilder();
        Document document = documentBuilder.parse(in);

        Element e = document.getDocumentElement();
        if (e.getTagName().equals("beast")) {

            concurrent = false;
            root = (XMLObject) convert(e, run);

        } else {
            throw new dr.xml.XMLParseException("Unknown root document element, " + e.getTagName());
        }

        return objectStore;
    }

    public XMLObject getRoot() {
        return root;
    }

    private Object convert(Element e, boolean run) throws XMLParseException {


        if (e.hasAttribute(IDREF)) {
            String idref = e.getAttribute(IDREF);
            XMLObject restoredXMLObject = (XMLObject) store.get(idref);
            if (restoredXMLObject == null) {
                throw new XMLParseException("Object with idref=" + idref + " has not been previously declared.");
            }
            if (verbose) System.out.println("  Restoring idref=" + idref);

            return new Reference(restoredXMLObject);
        } else {
            int repeats = 1;
            if (e.getTagName().equals(CONCURRENT)) {
                if (concurrent) throw new XMLParseException("Nested concurrent elements not allowed.");
                concurrent = true;

                threads = new Vector<Thread>();
            } else if (e.getTagName().equals("repeat")) {
                repeats = Integer.parseInt(e.getAttribute("count"));
            }

            XMLObject xo = new XMLObject(e /*, objectStore*/);

            String id = null;
            NodeList nodes = e.getChildNodes();
            for (int k = 0; k < repeats; k++) {
                for (int i = 0; i < nodes.getLength(); i++) {

                    Node child = nodes.item(i);
                    if (child instanceof Element) {

                        if (verbose) System.out.println("Parsing " + ((Element) child).getTagName());

                        xo.addChild(convert((Element) child, run));
                    } else if (child instanceof Text) {
                        // just add text as a child of type String object
                        String text = ((Text) child).getData().trim();
                        if (text.length() > 0) {
                            xo.addChild(text);
                        }
                    }
                }
            }
            if (e.hasAttribute(ID)) {
                id = e.getAttribute(ID);
            }

            if ((id != null) && store.get(id) != null) {
                throw new XMLParseException("Object with Id=" + id + " already exists");
            }

            XMLObjectParser parser = parserStore.get(xo.getName());

            Object obj = null;
            if (parser != null) {
                obj = parser.parseXMLObject(xo, id, objectStore);

                if (id != null && obj instanceof Identifiable) {
                    ((Identifiable) obj).setId(id);
                }

                xo.setNativeObject(obj);
            }

            if (id != null) {
                if (verbose) System.out.println("  Storing " + xo.getName() + " with id=" + id);

                store.put(id, xo);
            }

            if (run) {
                if (e.getTagName().equals(CONCURRENT)) {
                    for (int i = 0; i < xo.getChildCount(); i++) {
                        Object child = xo.getChild(i);
                        if (child instanceof Runnable) {
                            Thread thread = new Thread((Runnable) child);
                            thread.start();
                            threads.add(thread);
                        } else throw new XMLParseException("Concurrent element children must be runnable!");
                    }
                    concurrent = false;
                    // wait for all threads collected to die
                    for (Object thread1 : threads) {
                        waitForThread((Thread) thread1);
                    }
                } else if (obj instanceof Runnable && !concurrent) {
                    Thread thread = new Thread((Runnable) obj);
                    thread.start();
                    threads.add(thread);
                    waitForThread(thread);
                }
                threads.removeAllElements();
            }

            return xo;
        }
    }

    public class ArrayParser extends AbstractXMLObjectParser {

        public String getParserName() {
            return "array";
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            List<Object> list = new ArrayList<Object>();

            for (int i = 0; i < xo.getChildCount(); i++) {
                list.add(xo.getChild(i));
            }
            return list;
        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public String getParserDescription() {
            return "This element returns an array of the objects it contains.";
        }

        public Class getReturnType() {
            return Object[].class;
        }

        public XMLSyntaxRule[] getSyntaxRules() {

            return new XMLSyntaxRule[]{new ElementRule(Object.class, "Objects to be put in an array", 1, Integer.MAX_VALUE)};
        }

    }

    private void waitForThread(Thread thread) {
        // wait doggedly for thread to die
        while (thread.isAlive()) {
            try {
                thread.join();
            } catch (InterruptedException ie) {
                // DO NOTHING
            }
        }
    }

    //anonymous object store class
    private ObjectStore objectStore = new ObjectStore() {
        public Object getObjectById(Object uid) throws ObjectNotFoundException {
            XMLObject obj = (XMLObject) store.get(uid);
            if (obj == null) throw new ObjectNotFoundException("Object with uid=" + uid + " not found in ObjectStore");
            if (obj.hasNativeObject()) return obj.getNativeObject();
            return obj;
        }

        public boolean hasObjectId(Object uid) {
            Object obj = store.get(uid);
            return (obj != null);
        }

        public Set getIdSet() {
            return store.keySet();
        }

        public Collection getObjects() {
            return store.values();
        }

        public void addIdentifiableObject(Identifiable obj, boolean force) {

            String id = obj.getId();
            if (id == null) throw new IllegalArgumentException("Object must have a non-null identifier.");

            if (force) {
                store.put(id, obj);
            } else {
                if (store.get(id) == null) {
                    store.put(id, obj);
                }
            }
        }
    };

    private Hashtable<String, Object> store = new Hashtable<String, Object>();
    private TreeMap<String, XMLObjectParser> parserStore = new TreeMap<String, XMLObjectParser>(new ParserComparator());
    private boolean concurrent = false;
    private XMLObject root = null;

    private boolean verbose = false;

    public static class Utils {

        final static String GENERATIONS = "generations";
        final static String DAYS = "days";
        final static String MONTHS = "months";
        final static String YEARS = "years";
        final static String SUBSTITUTIONS = "substitutions";
        // Mutations has been replaced with substitutions...
        final static String MUTATIONS = "mutations";
        final static String UNKNOWN = "unknown";

        final static String UNITS = "units";

        /**
         * Converts a unit integer into a human-readable name.
         *
         * @param units
         * @return unit name
         */
        public static String getUnitString(Units.Type units) {
            switch (units) {
                case GENERATIONS:
                    return GENERATIONS;
                case DAYS:
                    return DAYS;
                case MONTHS:
                    return MONTHS;
                case YEARS:
                    return YEARS;
                    // Mutations has been replaced with substitutions...
                case SUBSTITUTIONS:
                    return SUBSTITUTIONS;
                default:
                    return UNKNOWN;
            }
        }


        /**
         * Throws a runtime exception if the element does not have
         * the given name.
         */
        public static void validateTagName(Element e, String name) throws XMLParseException {
            if (!e.getTagName().equals(name)) {
                throw new XMLParseException("Wrong tag name! Expected " + name + ", found " + e.getTagName() + ".");
            }
        }

        public static Units.Type getUnitsAttr(XMLObject xo) throws XMLParseException {

            Units.Type units = dr.evolution.util.Units.Type.GENERATIONS;
            if (xo.hasAttribute(UNITS)) {
                String unitsAttr = (String) xo.getAttribute(UNITS);
                if (unitsAttr.equals(YEARS)) {
                    units = dr.evolution.util.Units.Type.YEARS;
                } else if (unitsAttr.equals(MONTHS)) {
                    units = dr.evolution.util.Units.Type.MONTHS;
                } else if (unitsAttr.equals(DAYS)) {
                    units = dr.evolution.util.Units.Type.DAYS;
                } else if (unitsAttr.equals(SUBSTITUTIONS) || unitsAttr.equals(MUTATIONS)) {
                    units = dr.evolution.util.Units.Type.SUBSTITUTIONS;
                }
            }
            return units;
        }

        public static boolean hasAttribute(Element e, String name) {
            String attr = e.getAttribute(name);
            return ((attr != null) && !attr.equals(""));
        }

        /**
         * @return the first child element of the given name.
         */
        public static Element getFirstByName(Element parent, String name) {
            NodeList nodes = parent.getElementsByTagName(name);
            if (nodes.getLength() > 0) {
                return (Element) nodes.item(0);
            } else return null;
        }

    }

    class ParserComparator implements Comparator<String> {
        public int compare(String o1, String o2) {
            String name1 = o1.toUpperCase();
            String name2 = o2.toUpperCase();

            return name1.compareTo(name2);
        }
    }
}

   
