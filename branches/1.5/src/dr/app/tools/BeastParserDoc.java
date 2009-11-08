/*
 * BeastParserDoc.java
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

package dr.app.tools;

import dr.app.util.Arguments;
import dr.app.util.Utils;
import dr.app.beast.BeastParser;
import dr.app.beast.BeastVersion;
import dr.xml.*;
import dr.util.Version;

import java.io.File;
import java.io.PrintWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Properties;

public class BeastParserDoc {

    public final static String INDEX_HTML = "index.html";
    public final static String DEATAIL_HTML = "detail.html";

    private final static Version version = new BeastVersion();

    public final static String TITTLE = "BEAST " + version.getVersionString() + " Parser Library ("
            + version.getDateString() + ")";
    public final static String AUTHORS = "Alexei Drummond, Andrew Rambaut, Walter Xie";
    public final static String LINK1 = "http://beast.bio.ed.ac.uk/";
    public final static String LINK2 = "http://code.google.com/p/beast-mcmc/";

    private final BeastParser parser;

    public BeastParserDoc(BeastParser parser, String directory, boolean wikiFormat) throws java.io.IOException {
        this.parser = parser;

        File file = new File(directory);

        if (!file.exists()) {
            file.mkdir();
        }

        if (!file.isDirectory()) {
            throw new IllegalArgumentException(directory + " is not a directory!");
        }

        if (wikiFormat) {
            XMLDocumentationHandler handler = new WikiDocumentationHandler(parser);

            PrintWriter writer = new PrintWriter(new FileWriter(new File(directory, "xml_format.wiki")));

            System.out.println("Building element descriptions...");
            handler.outputElements(writer);
            System.out.println("done.");

            System.out.println("Building types table...");
            handler.outputIndex(writer);
            System.out.println("done.");

            writer.flush();
            writer.close();
        } else {
            XMLDocumentationHandler handler = new XMLDocumentationHandler(parser);

            System.out.println("Generate " + INDEX_HTML + " ...");
            PrintWriter writer = new PrintWriter(new FileWriter(new File(directory, INDEX_HTML)));

            handler.outputIndex(writer); // generate index.html
            System.out.println("done.");
            writer.flush();
            writer.close();

            System.out.println("Generate " + DEATAIL_HTML + " ...");
            writer = new PrintWriter(new FileWriter(new File(directory, DEATAIL_HTML)));

            handler.outputElements(writer);
            System.out.println("done.");

            writer.flush();
            writer.close();
        }
    }

//	private final void setup() throws XMLParseException {
    // add all the XMLObject parsers you need
//	}

    public static void printTitle() {

        System.out.println("+-----------------------------------------------\\");
        System.out.print("|");
        int n = 47 - TITTLE.length();
        int n1 = n / 2;
        int n2 = n1 + (n % 2);
        for (int i = 0; i < n1; i++) { System.out.print(" "); }
        System.out.print(TITTLE);
        for (int i = 0; i < n2; i++) { System.out.print(" "); }
        System.out.println("||");
        System.out.println("|   " + AUTHORS + " ||");
        System.out.println("|           " + LINK1 + "          ||");
        System.out.println("|      " + LINK2 + "     ||");
        System.out.println("\\-----------------------------------------------\\|");
        System.out.println(" \\-----------------------------------------------\\");
        System.out.println();
    }

    public static void printUsage(Arguments arguments) {

        arguments.printUsage("beastdoc", "<output-directory>");
        System.out.println();
        System.out.println("  Example: beastdoc ./doc");
        System.out.println();

    }

    //Main method
    public static void main(String[] args) throws java.io.IOException, XMLParseException {

        printTitle();

        Arguments arguments = new Arguments(
                new Arguments.Option[] {
                        new Arguments.Option("help", "option to print this message")
                });

        try {
            arguments.parseArguments(args);
        } catch (Arguments.ArgumentException ae) {
            System.out.println(ae);
            printUsage(arguments);
            return;
        }

        if (arguments.hasOption("help")) {
            printUsage(arguments);
            return;
        }

        String outputDirectory = null;

        String[] args2 = arguments.getLeftoverArguments();

        if (args2.length > 1) {
            System.err.println("Unknown option: " + args2[1]);
            System.err.println();
            printUsage(arguments);
            return;
        }

        if (args2.length > 0) {
            outputDirectory = args2[0];
        }

        if (outputDirectory == null) {
            // No input file name was given so throw up a dialog box...
            outputDirectory = Utils.getSaveFileName("BeastParserDoc " + version.getVersionString()
                    + " - Select output directory");
        }

        new BeastParserDoc(new BeastParser(new String[] {}, null, false, false), outputDirectory, false);

        System.exit(0);
    }
}

