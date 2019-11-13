/*
 * TwoParalogGeneConversionSubstitutionModel.java
 *
 * Copyright (c) 2002-2017 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.evomodelxml.branchmodel;

import dr.evolution.io.Importer;
import dr.evolution.io.NexusImporter;
import dr.evolution.tree.FlexibleTree;
import dr.evomodel.branchmodel.PairedParalogBranchModel;
import dr.evomodel.substmodel.SubstitutionModel;
import dr.evomodel.substmodel.geneconversion.PairedParalogGeneConversionSubstitutionModel;
import dr.inference.model.Parameter;
import dr.util.FileHelpers;
import dr.xml.*;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

/**
 * @author Xiang Ji
 * @author Jeff Thorne
 * @author Marc A. Suchard
 */
public class PairedParalogBranchModelParser extends AbstractXMLObjectParser {

    private final String NAME = "pairedParalogBranchModel";
    private final String SUBSTITUTION_NEXUS_FILE = "substAssignFile";

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        List<SubstitutionModel> substitutionModels = new ArrayList<>();

        for (PairedParalogGeneConversionSubstitutionModel geneConversionSubstitutionModel : xo.getAllChildren(PairedParalogGeneConversionSubstitutionModel.class)) {
            substitutionModels.add(geneConversionSubstitutionModel);
        }

        String nexusFileName = xo.getStringAttribute(SUBSTITUTION_NEXUS_FILE);
        final File file = FileHelpers.getFile(nexusFileName);
        NexusImporter importer;
        int[][] substitutionModelAssignment;
        FlexibleTree attributeTree;

        try {
            importer = new NexusImporter(new FileReader(file));
            attributeTree = (FlexibleTree) importer.importNextTree();
            attributeTree.resolveTree();
            substitutionModelAssignment = new int[attributeTree.getNodeCount()][];
            for (int i = 0; i < attributeTree.getNodeCount(); i++) {
                Object attribute = attributeTree.getNodeAttribute(attributeTree.getNode(i), "model");
                if (attribute instanceof Integer) {
                    substitutionModelAssignment[i] = new int[]{(int) attribute};
                } else if (attribute instanceof String) {
                    substitutionModelAssignment[i] = PairedParalogBranchModel.parseAssignmentString((String) attribute);
                } else {
                    throw new RuntimeException("Unrecognized assignment format");
                }
            }
            attributeTree.getRoot();
        } catch (FileNotFoundException e) {
            throw new XMLParseException(e.getMessage());
        } catch (Importer.ImportException e) {
            throw new XMLParseException(e.getMessage());
        } catch (IOException e) {
            throw new XMLParseException(e.getMessage());
        }


        if (substitutionModelAssignment.length != attributeTree.getNodeCount()) {
            throw new RuntimeException("Please specify the substitution model for every branch.");
        }

        List<Parameter> timeToDuplicaitonProportion = xo.getAllChildren(Parameter.class);

        return new PairedParalogBranchModel(NAME, substitutionModels, substitutionModelAssignment, timeToDuplicaitonProportion, attributeTree);
    }


    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(PairedParalogGeneConversionSubstitutionModel.class, 2, 2),
            new ElementRule(Parameter.class, "proportion of time until the duplication event on the duplication branch."),
            AttributeRule.newStringArrayRule(SUBSTITUTION_NEXUS_FILE),
    };

    @Override
    public String getParserDescription() {
        return null;
    }

    @Override
    public Class getReturnType() {
        return PairedParalogBranchModel.class;
    }

    @Override
    public String getParserName() {
        return NAME;
    }
}
