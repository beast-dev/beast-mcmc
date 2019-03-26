/*
 * BEASTGen.java
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

package dr.app.beastgen;

import dr.app.beauti.options.DateGuesser;
import dr.app.util.Arguments;
import dr.evolution.io.Importer;
import freemarker.template.*;

import java.io.*;
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class BEASTGen {

    public BEASTGen(DateGuesser guesser, Map argumentMap, String treeFileName, String taxonSetFileName ,String templateFileName, String inputFileName, String outputFileName) throws IOException {

        Configuration cfg = new Configuration();
        cfg.setObjectWrapper(new DefaultObjectWrapper());
        cfg.setDefaultEncoding("UTF-8");
        cfg.setTemplateExceptionHandler(TemplateExceptionHandler.DEBUG_HANDLER);
        cfg.setIncompatibleImprovements(new Version(2, 3, 20));  // FreeMarker 2.3.20

        Map root = null;
        try {
            root = constructDataModel(inputFileName, treeFileName, taxonSetFileName, guesser);
        } catch (Importer.ImportException ie) {
            System.err.println("Error importing file: " + ie.getMessage());
            System.exit(1);
        }

        root.putAll(argumentMap);

        Template template = cfg.getTemplate(templateFileName);

        try {
            Writer out = (outputFileName != null ?
                    new FileWriter(new File(outputFileName)) :
                    new OutputStreamWriter(System.out));

            template.process(root, out);
        } catch (TemplateException te) {
            System.err.println("Error processing template, " + templateFileName + ": " + te.getMessage());
            System.exit(1);
        }
    }

    private Map constructDataModel(String inputFileName, String treeFileName, String taxonSetFileName,DateGuesser guesser) throws IOException, Importer.ImportException {
        DataModelImporter importer = new DataModelImporter(guesser);

        Map root = importer.importFromFile(new File(inputFileName));

        if (treeFileName != null) {
            importer.importFromTreeFile(treeFileName, root);
        }
        if(taxonSetFileName !=null){
            importer.importTaxonSets(taxonSetFileName,root);
        }

        return root;
    }

    public static void centreLine(String line, int pageWidth) {
        int n = pageWidth - line.length();
        int n1 = n / 2;
        for (int i = 0; i < n1; i++) {
            System.out.print(" ");
        }
        System.out.println(line);
    }

    public static void printTitle() {
        System.out.println();
        centreLine("BEASTGen v1.0.1, 2013-2014", 60);
        centreLine("BEAST input file generator", 60);
        centreLine("Andrew Rambaut, University of Edinburgh", 60);
        System.out.println();
    }

    public static void printUsage(Arguments arguments) {

        arguments.printUsage("beastgen", "<template-file-name> <input-file-name> [<output-file-name>]");
        System.out.println();
        System.out.println("  Example: beastgen template.beast test.nex test.xml");
        System.out.println("  Example: beastgen -help");
        System.out.println();
    }

    public static void main(String[] args) {

        // There is a major issue with languages that use the comma as a decimal separator.
        // To ensure compatibility between programs in the package, enforce the US locale.
        Locale.setDefault(Locale.US);

        Arguments arguments = new Arguments(
                new Arguments.Option[]{
                        new Arguments.IntegerOption("date_order", "The order of the date field (negative numbers from last)"),
                        new Arguments.StringOption("date_prefix", "prefix", "A string that is the prefix to the date field"),
                        new Arguments.StringOption("date_regex", "regex", "A string that gives the regular expression to match the date"),
                        new Arguments.StringOption("date_format", "format", "A string that gives the date format for parsing"),
                        new Arguments.Option("date_precision", "Specifies the date is a variable precision yyyy-MM-dd format"),
                        new Arguments.StringOption("tree", "tree-file-name", "Read a tree from a file"),
                        new Arguments.StringOption("taxonSet", "taxonSet-file-name", "Read taxon sets from a tsv file with headers."),
                        new Arguments.StringOption("D", "\"key=value,key=value...\"", "Properties for exchange in templates"),
                        new Arguments.Option("version", "Print the version and credits and stop"),
                        new Arguments.Option("help", "Print this information and stop"),
                });


        try {
            arguments.parseArguments(args);
        } catch (Arguments.ArgumentException ae) {
            printTitle();
            System.out.println();
            System.out.println(ae.getMessage());
            System.out.println();
            printUsage(arguments);
            System.exit(1);
        }

        DateGuesser guesser = new DateGuesser();
//        guesser.guessDates = false;
//        guesser.guessType = DateGuesser.GuessType.ORDER;
//        guesser.fromLast = false;
//        guesser.order = 0;
//        guesser.prefix = null;
//        guesser.regex = null;
//        guesser.offset = 0.0;
//        guesser.unlessLessThan = 0.0;
//        guesser.offset2 = 0.0;
//        guesser.parseCalendarDates = false;
//        guesser.parseCalendarDatesAndPrecision = false;
//        guesser.calendarDateFormat = "yyyy-MM-dd";

        if (arguments.hasOption("date_order")) {
            guesser.guessDates = true;
            int order = arguments.getIntegerOption("date_order");
            if (order < 0) {
                guesser.order = 1 + order;
                guesser.fromLast = true;
            } else if (order > 0) {
                guesser.order = order - 1;
                guesser.fromLast = false;
            } else {
                guesser.order = 0;
                guesser.fromLast = false;
            }
        }

        if (arguments.hasOption("date_prefix")) {
            guesser.guessDates = true;
            guesser.prefix = arguments.getStringOption("date_prefix");
            guesser.guessType = DateGuesser.GuessType.PREFIX;
        }

        if (arguments.hasOption("date_regex")) {
            guesser.guessDates = true;
            guesser.regex = arguments.getStringOption("date_regex");
            guesser.guessType = DateGuesser.GuessType.REGEX;
        }

        if (arguments.hasOption("date_format")) {
            guesser.guessDates = true;
            guesser.calendarDateFormat = arguments.getStringOption("date_format");
            guesser.parseCalendarDates = true;
        }

        if (arguments.hasOption("date_precision")) {
            guesser.guessDates = true;
            guesser.parseCalendarDatesAndPrecision = true;
        }

        String treeFileName = null;
        if (arguments.hasOption("tree")) {
            treeFileName = arguments.getStringOption("tree");
        }

        String taxonSetFileName=null;
        if(arguments.hasOption("taxonSet")){
            taxonSetFileName=arguments.getStringOption("taxonSet");
        }


        Map argumentMap = new HashMap();

        if (arguments.hasOption("D")) {
            String properties = arguments.getStringOption("D");

            for (String property : properties.split("\\s*,\\s*")) {
                String[] keyValue = property.split("=");
                if (keyValue.length != 2) {
                    System.err.println("Properties should take the form: key=value");
                    System.exit(1);
                }
                String key = keyValue[0].trim();
                String value = keyValue[1].trim();
                if (key.isEmpty()) {
                    System.err.println("Properties should take the form: key=value");
                    System.exit(1);
                }
                argumentMap.put(key, value);
            }
        }

        if (arguments.hasOption("help")) {
            printTitle();
            printUsage(arguments);
            System.exit(0);
        }

        if (arguments.hasOption("version")) {
            printTitle();
        }


        String[] args2 = arguments.getLeftoverArguments();

        if (args2.length < 1) {
            printTitle();
            printUsage(arguments);
            System.exit(0);
        }

        if (args2.length < 2 || args2.length > 3) {
            System.err.println("Unknown option: " + args2[0]);
            System.err.println();
            printUsage(arguments);
            System.exit(1);
        }

        try {
            new BEASTGen(guesser, argumentMap, treeFileName, taxonSetFileName,args2[0], args2[1], (args2.length == 3 ? args2[2] : null));
        } catch (IOException e) {
            e.printStackTrace();
            System.exit(1);
        }

    }
}
