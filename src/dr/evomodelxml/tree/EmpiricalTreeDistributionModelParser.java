/*
 * EmpiricalTreeDistributionModelParser.java
 *
 * Copyright (c) 2002-2016 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.evomodelxml.tree;

import dr.util.FileHelpers;
import dr.xml.*;
import dr.evolution.io.NexusImporter;
import dr.evolution.io.Importer;
import dr.evolution.tree.Tree;
import dr.evolution.util.TaxonList;
import dr.evomodel.tree.EmpiricalTreeDistributionModel;

import java.io.*;
import java.util.logging.Logger;

/**
 * @author Andrew Rambaut
 *         <p/>
 *         Reads a list of trees from a NEXUS file.
 */
public class EmpiricalTreeDistributionModelParser extends AbstractXMLObjectParser {

    public static final String RATE_ATTRIBUTE_NAME = "rateAttribute";
    public static final String STARTING_TREE = "startingTree";
    public static final String ITERATE = "iterate";

    public String getParserName() {
        return EmpiricalTreeDistributionModel.EMPIRICAL_TREE_DISTRIBUTION_MODEL;
    }

    public String getParserDescription() {
        return "Read a list of trees from a NEXUS file.";
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        final String fileName = xo.getStringAttribute(FILE_NAME);

        int startingTree = xo.getAttribute(STARTING_TREE, -1); // default is random tree

        boolean iterate = xo.getAttribute(ITERATE, false); // default is random draw
        if (iterate && startingTree < 0) {
            startingTree = 0;
        }

        Logger.getLogger("dr.evomodel").info("Creating the empirical tree distribution model, '" + xo.getId() + "'");

        TaxonList taxa = (TaxonList)xo.getChild(TaxonList.class);

        final File file = FileHelpers.getFile(fileName);

        Tree[] trees = null;
        NexusImporter importer = null;
        try {
            FileReader reader = new FileReader(file);
            importer = new NexusImporter(reader);
            if (!iterate) {
                trees = importer.importTrees(taxa, true); // Re-order taxon numbers to original TaxonList order
                reader.close();
            }
        } catch (FileNotFoundException e) {
            throw new XMLParseException(e.getMessage());
        } catch (IOException e) {
            throw new XMLParseException(e.getMessage());
        } catch (Importer.ImportException e) {
            throw new XMLParseException(e.getMessage());
        }
        

        if (iterate) {
            Logger.getLogger("dr.evomodel").info("    Iterate over each tree from file, " + fileName);
            return new EmpiricalTreeDistributionModel(importer, startingTree);
        } else {
            Logger.getLogger("dr.evomodel").info("    Randomly jump between " + trees.length + " trees from file, " + fileName);
            return new EmpiricalTreeDistributionModel(trees, startingTree);
        }
    }

    public static final String FILE_NAME = "fileName";
//    public static final String BURNIN = "burnin";

    public XMLSyntaxRule[] getSyntaxRules() {
        return new XMLSyntaxRule[]{
                AttributeRule.newIntegerRule(STARTING_TREE, true),
                AttributeRule.newBooleanRule(ITERATE, true),
                new StringAttributeRule(FILE_NAME,
                        "The name of a NEXUS tree file"),
//                AttributeRule.newIntegerRule(BURNIN, true,
//                        "The number of trees to exclude"),
                new ElementRule(TaxonList.class),
        };
    }


    public Class getReturnType() {
        return EmpiricalTreeDistributionModel.class;
    }
}