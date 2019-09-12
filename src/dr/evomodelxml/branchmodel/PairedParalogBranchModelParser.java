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

import dr.evomodel.branchmodel.PairedParalogBranchModel;
import dr.evomodel.substmodel.SubstitutionModel;
import dr.evomodel.substmodel.geneconversion.DuplicationBranchSubstitutionModel;
import dr.evomodel.substmodel.geneconversion.PairedParalogGeneConversionSubstitutionModel;
import dr.evomodel.tree.TreeModel;
import dr.xml.*;

import java.util.ArrayList;
import java.util.List;

/**
 * @author Xiang Ji
 * @author Jeff Thorne
 * @author Marc A. Suchard
 */
public class PairedParalogBranchModelParser extends AbstractXMLObjectParser {

    private final String NAME = "pairedParalogBranchModel";
    private final String SUBST_ASSIGN = "substAssign";

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        int[] substitutionModelAssignment = xo.getIntegerArrayAttribute(SUBST_ASSIGN);
        List<SubstitutionModel> substitutionModels = new ArrayList<>();

        for (PairedParalogGeneConversionSubstitutionModel geneConversionSubstitutionModel : xo.getAllChildren(PairedParalogGeneConversionSubstitutionModel.class)) {
            substitutionModels.add(geneConversionSubstitutionModel);
        }

        for (DuplicationBranchSubstitutionModel duplicationBranchSubstitutionModel : xo.getAllChildren(DuplicationBranchSubstitutionModel.class)) {
            substitutionModels.add(duplicationBranchSubstitutionModel);
        }

        TreeModel tree = (TreeModel) xo.getChild(TreeModel.class);

        if (substitutionModelAssignment.length != tree.getNodeCount() - 1) {
            throw new RuntimeException("Please specify the substitution model for every branch.");
        }

        List<SubstitutionModel> assignedSubstitutionModels = new ArrayList<>();
        for (int i : substitutionModelAssignment) {
            assignedSubstitutionModels.add(substitutionModels.get(i));
        }

        return new PairedParalogBranchModel(NAME, assignedSubstitutionModels, tree);
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(PairedParalogGeneConversionSubstitutionModel.class, 1, Integer.MAX_VALUE),
            new ElementRule(DuplicationBranchSubstitutionModel.class, 1, Integer.MAX_VALUE),
            new ElementRule(TreeModel.class),
            AttributeRule.newIntegerArrayRule(SUBST_ASSIGN, false),
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
