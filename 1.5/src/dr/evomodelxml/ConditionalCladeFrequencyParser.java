/*
 * ConditionalCladeFrequencyParser.java
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

/**
 *
 */
package dr.evomodelxml;

import dr.evolution.io.Importer;
import dr.evolution.io.NewickImporter;
import dr.evolution.tree.Tree;
import dr.evomodel.tree.ConditionalCladeFrequency;
import dr.xml.*;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.Reader;

/**
 * @author Sebastian Hoehna
 *         <p/>
 *         This is the parser for the xml block in the BEAST input file for the conditional clade frequency statistic.
 *         It is used to estimate the posterior of a tree given a trace of trees.
 */
public class ConditionalCladeFrequencyParser extends AbstractXMLObjectParser {

    public final static String CONDITIONAL_CLADE_PROBABILITY_ANALYSIS = "ConditionalCladeProbabilityAnalysis";
    public final static String BURN_IN = "burnIn";
    public final static String EPSILON = "epsilon";
    public final static String MIN_CLADE_PROBABILITY = "minCladeProbability";
    public static final String FILE_NAME = "fileName";

    public final static String REFERENCE_TREE = "referenceTree";

    public String getParserName() {
        return CONDITIONAL_CLADE_PROBABILITY_ANALYSIS;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        try {
            Reader reader;

            // the tree file with the trace of trees. Usually an output of an previous BEAST run like *.trees
            String fileName = xo.getStringAttribute(FILE_NAME);
            String name;
            try {
                File file = new File(fileName);
                name = file.getName();
                String parent = file.getParent();

                if (!file.isAbsolute()) {
                    parent = System.getProperty("user.dir");
                }

                reader = new FileReader(new File(parent, name));
            } catch (FileNotFoundException fnfe) {
                throw new XMLParseException("File '" + fileName + "' can not be opened for " + getParserName() + " element.");
            }

            // the burn-in is used as the number of trees discarded
            int burnin = -1;
            if (xo.hasAttribute(BURN_IN)) {
                // leaving the burnin attribute off will result in 10% being used
                burnin = xo.getIntegerAttribute(BURN_IN);
            }

            // the epsilon value which represents the number of occurrences for every not observed clade
            double e = 1.0;
            if (xo.hasAttribute(EPSILON)) {
                // leaving the epsilon attribute off will result in 1.0 being used
                e = xo.getDoubleAttribute(EPSILON);
            }

            // not used yet
            double minCladeProbability = 0.5;
            if (xo.hasAttribute(MIN_CLADE_PROBABILITY)) {
                minCladeProbability = xo.getDoubleAttribute(MIN_CLADE_PROBABILITY);
            }

            // read the reference tree from a newick file
            Tree referenceTree = null;
            Reader refReader;
            if (xo.hasAttribute(REFERENCE_TREE)) {
                String referenceName = xo.getStringAttribute(REFERENCE_TREE);

                try {
                    File refFile = new File(referenceName);
                    String refName = refFile.getName();
                    String parent = refFile.getParent();

                    if (!refFile.isAbsolute()) {
                        parent = System.getProperty("user.dir");
                    }
                    refReader = new FileReader(new File(parent, refName));
                } catch (FileNotFoundException fnfe) {
                    throw new XMLParseException("File '" + fileName + "' can not be opened for " + getParserName() + " element.");
                }

                try {
                    NewickImporter importTree = new NewickImporter(refReader);
                    if (importTree.hasTree()) {
                        referenceTree = importTree.importNextTree();
                    }
                } catch (Importer.ImportException iee) {
                    throw new XMLParseException("Reference file '" + referenceName + "' is empty.");
                }
            }

            // creates a new analyzer object
            ConditionalCladeFrequency analysis = ConditionalCladeFrequency.analyzeLogFile(new Reader[]{reader}, e, burnin, true);

            // analyze the reference tree and prints its estimated posterior
            analysis.report(referenceTree);

            System.out.println();
            System.out.flush();

            return analysis;
        } catch (java.io.IOException ioe) {
            throw new XMLParseException(ioe.getMessage());
        }
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "Calculates posterior of a tree on a trace consisting of trees.";
    }

    public Class getReturnType() {
        return ConditionalCladeFrequency.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
            new StringAttributeRule(FILE_NAME, "name of a tree log file", "trees.log"),
            AttributeRule.newIntegerRule(BURN_IN, true, "The number of trees discarded because of the burn-in period."),
            AttributeRule.newDoubleRule(MIN_CLADE_PROBABILITY, true, "The frequency in % that a clade must achieve to be considered for the evaulation. Otherwise it's set to epsilon."),
            AttributeRule.newDoubleRule(EPSILON, true, "Epsilon is the default number of occurences for a clade if a clade wasn't observed. The default is 1.0"),
            AttributeRule.newStringRule(REFERENCE_TREE, false, "File in newick format containing a tree which is evaluated with the statistic.")
    };

}
