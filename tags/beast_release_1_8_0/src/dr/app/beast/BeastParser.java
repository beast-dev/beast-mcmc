/*
 * BeastParser.java
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

package dr.app.beast;

import dr.xml.PropertyParser;
import dr.xml.UserInput;
import dr.xml.XMLObjectParser;
import dr.xml.XMLParser;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.lang.reflect.Field;
import java.util.Iterator;
import java.util.List;
import java.util.Properties;

/**
 * @author Alexei Drummond
 * @author Andrew Rambaut
 * @author Walter Xie
 * @version $Id: BeastParser.java,v 1.76 2006/08/30 16:01:59 rambaut Exp $
 */
public class BeastParser extends XMLParser {

    public static final String RELEASE ="release";
    public static final String DEV = "development";
    public static final String PARSER_PROPERTIES_SUFFIX ="_parsers.properties";
    public String parsers;

    public BeastParser(String[] args, List<String> additionalParsers, boolean verbose, boolean parserWarnings, boolean strictXML) {
        super(parserWarnings, strictXML);

        setup(args);

        if (verbose) {
            System.out.println("Built-in parsers:");
            Iterator iterator = getParsers();
            while (iterator.hasNext()) {
                XMLObjectParser parser = (XMLObjectParser) iterator.next();
                System.out.println(parser.getParserName());
            }

        }

        // Try to find and load the additional 'core' parsers
        try {
            Properties properties = new Properties();
            properties.load(this.getClass().getResourceAsStream("beast.properties"));

            // get the parsers file prefix from the beast.properties file
            parsers = properties.getProperty("parsers");

            if (System.getProperty("parsers") != null) {
                // If a system property has been set then allow this to override the default
                // e.g. -Dparsers=development
                parsers = properties.getProperty("parsers");
            }

            if (parsers.equalsIgnoreCase(DEV)) {
                this.parserWarnings = true; // if dev, then auto turn on, otherwise default to turn off
            }

            // always load release_parsers.properties !!!
            loadProperties(this.getClass(), RELEASE + PARSER_PROPERTIES_SUFFIX, verbose, this.parserWarnings, false);

            // suppose to load developement_parsers.properties
            if (parsers != null && (!parsers.equalsIgnoreCase(RELEASE))) {
                // load the development parsers
                if (parsers.equalsIgnoreCase(DEV)) {
                    System.out.println("\nLoading additional development parsers from " + parsers + PARSER_PROPERTIES_SUFFIX
                            + ", which is additional set of parsers only available for development version ...");
                }
                loadProperties(this.getClass(), parsers + PARSER_PROPERTIES_SUFFIX, verbose, this.parserWarnings, true);
            }
            // load additional parsers
            if (additionalParsers != null) {
                for (String addParsers : additionalParsers) {
                    loadProperties(this.getClass(), addParsers + PARSER_PROPERTIES_SUFFIX, verbose, verbose, true);
                }
            }
        } catch (IOException e) {
            e.printStackTrace();
        }

        // Now search the package hierarchy for 'beast.properties' files.
//        try {
//            loadProperties(this.getClass(), verbose);
//        } catch (IOException e) {
//            e.printStackTrace();
//        }
    }

    /**
     * Load the parser for *.properties file
     * @param c               BeastParser
     * @param parsersFile     parser file name, (*.properties)
     * @param verbose         verbose
     * @param parserWarning   parserWarning
     * @param canReplace      can this new loaded parser to replace old one with the same name
     * @throws IOException    IOException
     */
    private void loadProperties(Class c, String parsersFile, boolean verbose, boolean parserWarning, boolean canReplace) throws IOException {

        if (verbose) {
            if (parsersFile.equalsIgnoreCase(RELEASE + PARSER_PROPERTIES_SUFFIX)) {
                System.out.println("\nAlways loading " + parsersFile + ":");
            } else {
                System.out.println("\n\nLoading additional parsers (" + parsersFile + "):");
            }
        }
        final InputStream stream = c.getResourceAsStream(parsersFile);
        if (stream == null) {
            throw new RuntimeException("Parsers file not found: " + parsersFile);
        }
        BufferedReader reader = new BufferedReader(new InputStreamReader(stream));
        String line = reader.readLine();

        while (line != null) {
            if (verbose && line.trim().startsWith("#")) System.out.println(line);

            if (line.trim().length() > 0 && !line.trim().startsWith("#")) {
                try {
                    if (line.contains("Vector")) {
                        System.out.println("");
                    }
                    Class parser = Class.forName(line);
                    if (XMLObjectParser.class.isAssignableFrom(parser)) {
                        // if this class is an XMLObjectParser then create an instance
                        boolean replaced = addXMLObjectParser((XMLObjectParser) parser.newInstance(), canReplace);
                        if (verbose) {
                            System.out.println((replaced ? "Replaced" : "Loaded") + " parser: " + parser.getName());
                        } else if (parserWarning && replaced) {
                            System.out.println("WARNING: parser - " + parser.getName() + " in " + parsersFile +" is duplicated, "
                                    + "which is REPLACING the same parser loaded previously.\n");
                        }
                    } else {
                        boolean parserFound = false;
                        // otherwise look for a static member which is an instance of XMLObjectParser
                        Field[] fields = parser.getDeclaredFields();
                        for (Field field : fields) {
                            if (XMLObjectParser.class.isAssignableFrom(field.getType())) {
                                try {
                                    boolean replaced = addXMLObjectParser((XMLObjectParser) field.get(null), canReplace);
                                    if (verbose) {
                                        System.out.println((replaced ? "Replaced" : "Loaded") + " parser: "
                                                + parser.getName() + "." + field.getName());
                                    } else if (parserWarning && replaced) {
                                        System.out.println("WARNING: parser - " + parser.getName() + " in " + parsersFile +" is duplicated, "
                                                + "which is REPLACING the same parser loaded previously.\n");
                                    }
                                } catch (IllegalArgumentException iae) {
                                    System.err.println("Failed to install parser: " + iae.getMessage());
                                }
                                parserFound = true;
                            }
                        }

                        if (!parserFound) {
                            throw new IllegalArgumentException(parser.getName() + " is not of type XMLObjectParser " +
                                    "and doesn't contain any static members of this type");
                        }
                    }

                } catch (Exception e) {
                    System.err.println("\nFailed to load parser: " + e.getMessage());
                    System.err.println("line = " + line + "\n");
                }
            }
            line = reader.readLine();
        }

        if (verbose) {
            System.out.println("load " + parsersFile + " successfully.\n");
        }
    }

    private void setup(String[] args) {

        for (int i = 0; i < args.length; i++) {
            storeObject(Integer.toString(i), args[i]);
        }

        // built-in parsers

        addXMLObjectParser(new PropertyParser());
        addXMLObjectParser(UserInput.STRING_PARSER);
        addXMLObjectParser(UserInput.DOUBLE_PARSER);
        addXMLObjectParser(UserInput.INTEGER_PARSER);

        addXMLObjectParser(new dr.xml.AttributeParser());
        addXMLObjectParser(new dr.xml.AttributesParser());

        addXMLObjectParser(new dr.inference.model.StatisticParser());
        addXMLObjectParser(new dr.inference.model.ParameterParser());

        //**************** all other parsers are read at runtime from property lists *********************
    }
}

