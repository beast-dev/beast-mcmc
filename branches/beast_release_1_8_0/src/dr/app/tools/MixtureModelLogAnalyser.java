/*
 * MixtureModelLogAnalyser.java
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
 *
 */

package dr.app.tools;

import dr.app.beast.BeastVersion;
import dr.app.util.Arguments;
import dr.app.util.Utils;
import dr.inference.trace.PathSamplingAnalysis;
import dr.inference.trace.TraceException;
import dr.util.Attribute;
import dr.util.FileHelpers;
import dr.util.Version;
import dr.xml.AbstractXMLObjectParser;
import dr.xml.AttributeRule;
import dr.xml.ElementRule;
import dr.xml.StringAttributeRule;
import dr.xml.XMLObject;
import dr.xml.XMLObjectParser;
import dr.xml.XMLParseException;
import dr.xml.XMLSyntaxRule;

import java.io.*;
import java.util.*;

/**
 * @author Wai Lok Sibon Li
 */

public class MixtureModelLogAnalyser {

    private final static Version version = new BeastVersion();
    
    public static final String MIXTURE_MODEL_LOG_ANALYSER = "mixtureModelLogAnalyser";
    public static final String BURNIN = "burnin";
    public static final String DISCRETE_VARIABLE = "discreteVariable";

    public MixtureModelLogAnalyser(int burnin, String inputFileName, String outputFileName, String discreteVariableName
        ) throws IOException, TraceException {

        File parentFile = new File(inputFileName);

        if (parentFile.isDirectory()) {
            System.out.println("Analysing all log files below directory: " + inputFileName);
        } else if (parentFile.isFile()) {
            System.out.println("Analysing log file: " + inputFileName);
        } else {
            System.err.println(inputFileName + " does not exist!");
            System.exit(0);
        }

        if (outputFileName != null) {
            FileOutputStream outputStream = new FileOutputStream(outputFileName);
            System.setOut(new PrintStream(outputStream));
        }

        analyze(parentFile, burnin, discreteVariableName);
    }

    /**
     * Recursively analyzes log files.
     *
     * @param file       the file to analyze (if this is a directory then the files within it are analyzed)
     * @param burnin     the burnin to use
     * @param discreteVariableName  tag for the name of the discrete variable
     * @throws dr.inference.trace.TraceException
     *          if the trace file is in the wrong format or corrupted
     */
    private void analyze(File file, int burnin, String discreteVariableName) throws TraceException {

        if (file.isFile()) {
            try {

                String name = file.getCanonicalPath();

                report(name, burnin, discreteVariableName);

            } catch (IOException e) {
                //e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
            }
        } else {
            File[] files = file.listFiles();
            for (File f : files) {
                if (f.isDirectory()) {
                    analyze(f, burnin, discreteVariableName);
                } else if (f.getName().endsWith(".log")) {
                    analyze(f, burnin, discreteVariableName);
                }
            }
        }
    }


    /**
     * Recursively analyzes trees files.
     *
     * @param name                  the file to analyze (if this is a directory then the files within it are analyzed)
     * @param burnin                the burnin to use
     * @param discreteVariableName  tag for the name of the discrete variable
     */
    private void report(String name, int burnin, String discreteVariableName) {

        try {
            FileReader fileReader = new FileReader(new File(name));
            BufferedReader br = new BufferedReader(fileReader);

            String line;
            String previousLine =null;
            int discreteVariableIndex = -1;
            Hashtable<String, Integer> hash = new Hashtable<String, Integer>();
            while((line=br.readLine()) != null) {

                if(line.matches("\\d+\\t.+")) {
                    String[] split = line.split("\t");
                    if(split[0].equals("0")) {
                        String[] headerSplit = previousLine.split("\t");
                        for(int i=0; i<headerSplit.length; i++) {
                            String s = headerSplit[i];
                            if(s.equals(discreteVariableName)) {
                                discreteVariableIndex = i;
                            }
                        }
                    }

                    if(Integer.parseInt(split[0])>=burnin) {
                        if(hash.containsKey(split[discreteVariableIndex])) {
                            hash.put(split[discreteVariableIndex], new Integer(hash.get(split[discreteVariableIndex]) + 1));
                        }
                        else {
                            hash.put(split[discreteVariableIndex], new Integer(1));
                        }
                    }
                }

                previousLine = line;
            }

            Vector<String> v = new Vector(hash.keySet());
            Collections.sort(v);


            System.out.print("name" + "\t");
            for(String s : v) {

                double state = Double.parseDouble(s);

                if(Math.floor(state)!=state) {
                    throw new RuntimeException("State is not a whole number: " + state);
                }

                System.out.print("state " + (int) state + "\t");
            }
            System.out.print("\n" + name + "\t");
            for (String aV : v) {
                String element = aV;
                System.out.print(hash.get(element) + "\t");
            }
            System.out.println();
            br.close();
        } catch (IOException e) {
            System.err.println("Error Parsing Input log: " + e.getMessage());
        }
    }
    
    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return MIXTURE_MODEL_LOG_ANALYSER;
        }
        
        public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        	
        	try {
        		
        		String inputFileName = xo.getStringAttribute(FileHelpers.FILE_NAME);
        		String discreteVariableName = xo.getStringAttribute(DISCRETE_VARIABLE);
        	
        		int burninLength = 0;
        		if (xo.hasAttribute(BURNIN)) {
        			burninLength = xo.getIntegerAttribute(BURNIN);
        		}
            
        		MixtureModelLogAnalyser mixtureModel = new MixtureModelLogAnalyser(burninLength, inputFileName, null, discreteVariableName);
            
        		return mixtureModel;
            
        	} catch (IOException ioe) {
        		throw new XMLParseException(ioe.getMessage());
       		} catch (TraceException te) {
       			throw new XMLParseException(te.getMessage());
       		}

        }
        
      //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public String getParserDescription() {
            return "Performs posterior probabilities calculations.";
        }

        public Class getReturnType() {
            return MixtureModelLogAnalyser.class;
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private final XMLSyntaxRule[] rules = {
                new StringAttributeRule(FileHelpers.FILE_NAME,
                        "The traceName of a BEAST log file (can not include trees, which should be logged separately"),
                new StringAttributeRule(DISCRETE_VARIABLE,
                        "The tag name of the discrete variable"),       
                AttributeRule.newIntegerRule(BURNIN, true),
        };
        
    };

    public static void printTitle() {
        System.out.println();
        centreLine("MixtureModelLogAnalyser " + version.getVersionString() + ", " + version.getDateString(), 60);
        centreLine("MCMC Output analysis", 60);
        centreLine("by", 60);
        centreLine("Andrew Rambaut and Alexei J. Drummond", 60);
        System.out.println();
        centreLine("Institute of Evolutionary Biology", 60);
        centreLine("University of Edinburgh", 60);
        centreLine("a.rambaut@ed.ac.uk", 60);
        System.out.println();
        centreLine("Department of Computer Science", 60);
        centreLine("University of Auckland", 60);
        centreLine("alexei@cs.auckland.ac.nz", 60);
        System.out.println();
        System.out.println();
    }

    public static void centreLine(String line, int pageWidth) {
        int n = pageWidth - line.length();
        int n1 = n / 2;
        for (int i = 0; i < n1; i++) {
            System.out.print(" ");
        }
        System.out.println(line);
    }


    public static void printUsage(Arguments arguments) {

        arguments.printUsage("mixturemodelloganalyser", "[-burnin <burnin>][<input-file-name> [<output-file-name>]]");
        System.out.println();
        System.out.println("  Example: treelengthfinder test.log");
        System.out.println("  Example: treelengthfinder -burnin 10000 -discreteVariable branchRates.distributionIndex trees.log out.txt");
        System.out.println();

    }

    //Main method
    public static void main(String[] args) throws IOException, TraceException {

        printTitle();

        Arguments arguments = new Arguments(
                new Arguments.Option[]{
                        new Arguments.IntegerOption("burnin", "the number of states to be considered as 'burn-in'"),
                        new Arguments.StringOption("discreteVariable", "variable_name", "indicates the name of a variable that is actually discrete in nature"),
                        new Arguments.Option("help", "option to print this message")
                });

        try {
            arguments.parseArguments(args);
        } catch (Arguments.ArgumentException ae) {
            System.out.println(ae);
            printUsage(arguments);
            System.exit(1);
        }

        if (arguments.hasOption("help")) {
            printUsage(arguments);
            System.exit(0);
        }

        int burnin = -1;
        if (arguments.hasOption("burnin")) {
            burnin = arguments.getIntegerOption("burnin");
        }

        String discreteVariableName = null; 
        if (arguments.hasOption("discreteVariable")) {
            discreteVariableName = arguments.getStringOption("discreteVariable");
        }

        String inputFileName = null;
        String outputFileName = null;

        String[] args2 = arguments.getLeftoverArguments();

        if (args2.length > 2) {
            System.err.println("Unknown option: " + args2[2]);
            System.err.println();
            printUsage(arguments);
            System.exit(1);
        }

        if (args2.length > 0) {
            inputFileName = args2[0];
        }
        if (args2.length > 1) {
            outputFileName = args2[1];
        }

        if (inputFileName == null) {
            // No input file name was given so throw up a dialog box...
            inputFileName = Utils.getLoadFileName("MixtureModelLogAnalyser " + version.getVersionString() + " - Select log file to analyse");
        }

        if(burnin==-1) {
            System.out.println("Enter number of trees to burn-in (integer): ");
            BufferedReader br = new BufferedReader(new InputStreamReader(System.in));
            burnin = Integer.parseInt(br.readLine());
        }
        if(discreteVariableName==null) {
            System.out.println("Enter name of discrete variable: ");
            BufferedReader br = new BufferedReader(new InputStreamReader(System.in));
            discreteVariableName = br.readLine();
        }

        new MixtureModelLogAnalyser(burnin, inputFileName, outputFileName, discreteVariableName/*, !shortReport, hpds, ess, stdErr, marginalLikelihood*/);

        System.exit(0);
    }
}