/*
 * XMLParser.java
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

import dr.inference.model.Likelihood;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inferencexml.loggers.LoggerParser;
import dr.util.*;
import org.w3c.dom.*;
import org.xml.sax.InputSource;
import org.xml.sax.SAXException;
import org.xml.sax.SAXParseException;
import org.xml.sax.helpers.DefaultHandler;

import java.io.*;
import java.util.*;

public class XMLParser {

    public static final String ID = XMLObject.ID;
    public static final String IDREF = "idref";
    public static final String CONCURRENT = "concurrent";
    public static final String VERSION = "version";

    private Vector<Thread> threads = new Vector<Thread>();
    protected boolean strictXML;
    protected boolean parserWarnings;

    // The software version
    private final Version version;

    public XMLParser(boolean verbose, boolean parserWarnings, boolean strictXML, Version version) {
        this.verbose = verbose;
        this.parserWarnings = parserWarnings;
        this.strictXML = strictXML;
        addXMLObjectParser(new ArrayParser(), false);
        addXMLObjectParser(Report.PARSER, false);

        this.version = version;
    }

    public void addXMLObjectParser(XMLObjectParser parser) {
        addXMLObjectParser(parser, false);
    }

    public boolean addXMLObjectParser(XMLObjectParser parser, boolean canReplace) {

        boolean replaced = false;
        String[] parserNames = parser.getParserNames();

        for (String parserName : parserNames) {
            XMLObjectParser oldParser = parserStore.get(parserName);
            if (oldParser != null) {
                if (!canReplace) {
                    throw new IllegalArgumentException("New parser (" + parser.getParserName()
                            + ") in {" + parser.getReturnType() + "} cannot replace existing parser ("
                            + oldParser.getParserName() + ") in {" + oldParser.getReturnType() + "}");
                } else {
                    replaced = true;
                }
            }
            parserStore.put(parserName, parser);
        }

        return replaced;
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

        XMLObject xo = new XMLObject(null, null /*, objectStore*/);
        xo.setNativeObject(object);
        objectStore.put(name, xo);
    }

    /**
     * An alternative parser that parses until it finds an object of the given
     * class and then returns it.
     *
     * @param reader the reader
     * @param target the target class
     * @return
     * @throws java.io.IOException
     * @throws org.xml.sax.SAXException
     * @throws dr.xml.XMLParseException
     * @throws javax.xml.parsers.ParserConfigurationException
     *
     */
    public Object parse(Reader reader, Class target)
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
            return convert(e, target, null, false, true);

        } else {
            throw new dr.xml.XMLParseException("Unknown root document element, " + e.getTagName());
        }
    }

    public Map<String, XMLObject> parse(Reader reader, boolean run)
            throws java.io.IOException,
            org.xml.sax.SAXException,
            dr.xml.XMLParseException,
            javax.xml.parsers.ParserConfigurationException {

        InputSource in = new InputSource(reader);
        javax.xml.parsers.DocumentBuilderFactory documentBuilderFactory = javax.xml.parsers.DocumentBuilderFactory.newInstance();

        javax.xml.parsers.DocumentBuilder documentBuilder = documentBuilderFactory.newDocumentBuilder();
        documentBuilder.setErrorHandler(new MyErrorHandler());
        Document document = documentBuilder.parse(in);

        Element e = document.getDocumentElement();
        if (e.getTagName().equals("beast")) {
            // If the 'version' is attribute is present then check it is not an more recent version...
            if (e.hasAttribute(VERSION)) {
                String xmlVersion = e.getAttribute(VERSION);
                if (version != null && Version.Utils.isMoreRecent(xmlVersion, version.getVersion())) {
                   throw new XMLParseException("The version of BEAUti that generated this XML (" + xmlVersion + ") is more recent than the version of BEAST running it (" + version.getVersion() + "). This may be incompatible and cause unpredictable errors.");
                }
            }

            concurrent = false;
            root = (XMLObject) convert(e, null, null, run, true);

        } else {
            throw new dr.xml.XMLParseException("Unknown root document element, " + e.getTagName());
        }

        return objectStore;
    }

    private class MyErrorHandler extends DefaultHandler {
        public void warning(SAXParseException e) throws SAXException {
            System.out.println("Warning: ");
            printInfo(e);
        }

        public void error(SAXParseException e) throws SAXException {
            System.out.println("Error: ");
            printInfo(e);
        }

        public void fatalError(SAXParseException e) throws SAXException {
            System.out.println("Fatal error: ");
            printInfo(e);
        }

        private void printInfo(SAXParseException e) {
//            System.out.println("\tPublic ID: " + e.getPublicId());
//            System.out.println("\tSystem ID: " + e.getSystemId());
            System.out.println("\tLine number: " + e.getLineNumber());
            System.out.println("\tColumn number: " + e.getColumnNumber());
            System.out.println("\tError message: " + e.getMessage());
        }
    }

    public XMLObject getRoot() {
        return root;
    }

    private Object convert(Element e, Class target, XMLObject parent, boolean run, boolean doParse) throws XMLParseException {

        int index = -1;

        if (e.hasAttribute(IDREF)) {

            String idref = e.getAttribute(IDREF);

            if (e.hasAttribute("index")) {
                index = Integer.parseInt(e.getAttribute("index"));
            }
            if ((e.getAttributes().getLength() > 1 || e.getChildNodes().getLength() > 1) && index == -1) {
                throw new XMLParseException("Object with idref=" + idref + " must not have other content or attributes (or perhaps it was not intended to be a reference?).");
            }


            XMLObject restoredXMLObject = objectStore.get(idref);
            if (index != -1) {

                if (restoredXMLObject.getNativeObject() instanceof List) {

                    restoredXMLObject = new XMLObject(restoredXMLObject, index);
                } else {
                    throw new XMLParseException("Trying to get indexed object from non-list");
                }
            }

            if (restoredXMLObject == null) {
                throw new XMLParseException("Object with idref=" + idref + " has not been previously declared.");
            }

            if (restoredXMLObject.getNativeObject() == null) {
                throw new XMLParseException("Object with idref=" + idref + " has not been parsed.");
            }

            XMLObjectParser parser = parserStore.get(e.getTagName());
            boolean classMatch = parser != null && parser.getReturnType().isAssignableFrom(restoredXMLObject.getNativeObject().getClass());

            if (!e.getTagName().equals(restoredXMLObject.getName()) && !classMatch) {
                String msg = "Element named " + e.getTagName() + " with idref=" + idref +
                        " does not match stored object with same id and tag name " + restoredXMLObject.getName();
                if (strictXML) {
                    throw new XMLParseException(msg);
                } else if (parserWarnings) {
//                    System.err.println("WARNING: " + msg);
                    java.util.logging.Logger.getLogger("dr.xml").warning(msg);
                }
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

            XMLObject xo = new XMLObject(e, parent);

            final XMLObjectParser parser = doParse ? parserStore.get(xo.getName()) : null;

            String id = null;
            NodeList nodes = e.getChildNodes();
            for (int k = 0; k < repeats; k++) {
                for (int i = 0; i < nodes.getLength(); i++) {

                    Node child = nodes.item(i);
                    if (child instanceof Element) {

                        final Element element = (Element) child;
                        final String tag = element.getTagName();
                        if (verbose) System.out.println("Parsing " + tag);

                        // don't parse elements that may be legal here with global parsers
                        final boolean parseIt = parser == null || !parser.isAllowed(tag);
                        Object xoc = convert(element, target, xo, run, parseIt);
                        xo.addChild(xoc);

                        if (target != null && xoc instanceof XMLObject) {
                            Object obj = ((XMLObject) xoc).getNativeObject();
                            if (obj != null && target.isInstance(obj)) {
                                return obj;
                            }
                        }

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

            if ((id != null) && objectStore.get(id) != null) {
                throw new XMLParseException("Object with Id=" + id + " already exists");
            }

            Object obj = null;
            if (parser != null) {
                obj = parser.parseXMLObject(xo, id, objectStore, strictXML);

                if (obj instanceof Identifiable) {
                    ((Identifiable) obj).setId(id);
                }

                if (obj instanceof Citable) {
                    addCitable((Citable)obj);
                }

                if (obj instanceof Likelihood) {
                    Likelihood.FULL_LIKELIHOOD_SET.add((Likelihood) obj);
                } else if (obj instanceof Model) {
                    Model.FULL_MODEL_SET.add((Model) obj);
                } else if (obj instanceof Parameter) {
                    Parameter.FULL_PARAMETER_SET.add((Parameter) obj);
                }

                xo.setNativeObject(obj);
            }

            if (id != null) {
                if (verbose) System.out.println("  Storing " + xo.getName() + " with id=" + id);

                objectStore.put(id, xo);
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

                    executingRunnable();

                    if (obj instanceof Spawnable && !((Spawnable) obj).getSpawnable()) {
                        ((Spawnable) obj).run();
                    } else {
                        Thread thread = new Thread((Runnable) obj);
                        thread.start();
                        threads.add(thread);
                        waitForThread(thread);
                    }
                }
                threads.removeAllElements();
            }

            return xo;
        }
    }

    protected void executingRunnable() {
        // do nothing - for overriding by subclasses
    }

    public Map<Pair<String, String>, List<Citation>> getCitationStore() {
        return citationStore;
    }

    public static FileReader getFileReader(XMLObject xo, String attributeName) throws XMLParseException {
        if (xo.hasAttribute(attributeName)) {
            final File inFile = getFileHandle(xo, attributeName);
            try {
                return new FileReader(inFile);
            } catch (FileNotFoundException e) {
                throw new XMLParseException("Input file " + inFile.getName()
                        + " was not found in the working directory");
            }
        }
        throw new XMLParseException("Error reading input file in " + xo.getId());
    }


    /**
     * Get filename and path from BEAST XML object
     *
     * @param xo
     * @return
     */
    private static File getFileHandle(XMLObject xo, String attributeName) throws XMLParseException {
        String fileName = xo.getStringAttribute(attributeName);

        // Check to see if a filename prefix has been specified, check it doesn't contain directory
        // separator characters and then prefix it.
        final String fileNamePrefix = System.getProperty("file.name.prefix");
        final String fileSeparator = System.getProperty("file.separator");
        if (fileNamePrefix != null) {
            if (fileNamePrefix.trim().length() == 0 || fileNamePrefix.contains(fileSeparator)) {
                throw new XMLParseException("The specified file name prefix is illegal.");
            }
        }

        final String fileRankPostfix = System.getProperty("mpi.rank.postfix");
        if (fileRankPostfix != null) {
            if (fileName.endsWith(".log")) {
                fileName = fileName.substring(0, fileName.length() - 4) + fileRankPostfix + ".log";
            }
            if (fileName.endsWith(".trees")) {
                fileName = fileName.substring(0, fileName.length() - 6) + fileRankPostfix + ".trees";
            }
        }

        return FileHelpers.getFile(fileName, fileNamePrefix);
    }

    /**
     * Allow a file relative to beast xml file with a prefix of ./
     *
     * @param xo         element
     * @param parserName for error messages
     * @return Print writer from fileName attribute in the given XMLObject
     * @throws XMLParseException if file can't be created for some reason
     */

    public static PrintWriter getFilePrintWriter(XMLObject xo, String parserName) throws XMLParseException {
        return getFilePrintWriter(xo, parserName, FileHelpers.FILE_NAME);
    }

    public static PrintWriter getFilePrintWriter(XMLObject xo, String parserName, String attributeName) throws XMLParseException {

        if (xo.hasAttribute(attributeName)) {
            File logFile = getLogFile(xo, attributeName);

            try {
                return new PrintWriter(new FileOutputStream(logFile));
            } catch (FileNotFoundException fnfe) {
                throw new XMLParseException("File '" + logFile.getAbsolutePath() +
                        "' can not be opened for " + parserName + " element.");
            }

        }
        return new PrintWriter(System.out);
    }

    public static File getLogFile(XMLObject xo, String attributeName) throws XMLParseException {
        final File logFile = getFileHandle(xo, attributeName);
        boolean allowOverwrite = false;

        if (xo.hasAttribute(LoggerParser.ALLOW_OVERWRITE_LOG)) {
            allowOverwrite = xo.getBooleanAttribute(LoggerParser.ALLOW_OVERWRITE_LOG);
        }

        // override with a runtime set System Property
        if (System.getProperty("log.allow.overwrite") != null) {
            allowOverwrite = Boolean.parseBoolean(System.getProperty("log.allow.overwrite", "false"));
        }

        if (logFile.exists() && !allowOverwrite) {
            throw new XMLParseException("\nThe log file " + logFile.getName() + " already exists in the working directory." +
                    "\nTo allow it to be overwritten, use the '-overwrite' command line option when running" +
                    "\nBEAST or select the option in the Run Options dialog box as appropriate.");
        }


        return logFile;
    }

    public Map<String, XMLObject> getObjectStore() {
        return objectStore;
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

//    //anonymous object store class
//    private final ObjectStore objectStore = new ObjectStore() {
//        public Object getObjectById(Object uid) throws ObjectNotFoundException {
//            XMLObject obj = (XMLObject) store.get(uid);
//            if (obj == null) throw new ObjectNotFoundException("Object with uid=" + uid + " not found in ObjectStore");
//            if (obj.hasNativeObject()) return obj.getNativeObject();
//            return obj;
//        }
//
//        public boolean hasObjectId(Object uid) {
//            Object obj = store.get(uid);
//            return (obj != null);
//        }
//
//        public Set getIdSet() {
//            return store.keySet();
//        }
//
//        public Collection getObjects() {
//            return store.values();
//        }
//
////        public void addIdentifiableObject(Identifiable obj, boolean force) {
////
////            String id = obj.getId();
////            if (id == null) throw new IllegalArgumentException("Object must have a non-null identifier.");
////
////            if (force) {
////                store.put(id, obj);
////            } else {
////                if (store.get(id) == null) {
////                    store.put(id, obj);
////                }
////            }
////        }
//    };

    public void addCitable(Citable citable) {
        // remove 'In prep' citations
        List<Citation> citationList = new LinkedList<Citation>();
        for (Citation citation : citable.getCitations()) {
            if (citation.getStatus() != Citation.Status.IN_PREPARATION) {
                citationList.add(citation);
            }
        }
        if (citationList.size() > 0) {
            Pair<String, String> pair = new Pair<String, String>(citable.getCategory().toString(),
                    citable.getDescription());
            citationStore.put(pair, citationList);
        }
    }

    //    private final Hashtable<String, XMLObject> store = new Hashtable<String, XMLObject>();
    private final Map<String, XMLObjectParser> parserStore = new TreeMap<String, XMLObjectParser>(new ParserComparator());
    private final Map<String, XMLObject> objectStore = new LinkedHashMap<String, XMLObject>();
    private final Map<Pair<String, String>, List<Citation>> citationStore = new LinkedHashMap<Pair<String, String>, List<Citation>>();
    private boolean concurrent = false;
    private XMLObject root = null;

    private boolean verbose = false;

    public static class Utils {

        /**
         * Throws a runtime exception if the element does not have
         * the given name.
         */
        public static void validateTagName(Element e, String name) throws XMLParseException {
            if (!e.getTagName().equals(name)) {
                throw new XMLParseException("Wrong tag name! Expected " + name + ", found " + e.getTagName() + ".");
            }
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


