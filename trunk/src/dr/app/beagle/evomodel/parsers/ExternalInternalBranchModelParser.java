/*
 * BranchSpecificBranchModelParser.java
 *
 * Copyright (C) 2002-2012 Alexei Drummond, Andrew Rambaut & Marc A. Suchard
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

package dr.app.beagle.evomodel.parsers;

import dr.app.beagle.evomodel.branchmodel.BranchSpecificBranchModel;
import dr.app.beagle.evomodel.branchmodel.ExternalInternalBranchModel;
import dr.app.beagle.evomodel.substmodel.SubstitutionModel;
import dr.evolution.tree.Tree;
import dr.evolution.util.Taxa;
import dr.evolution.util.TaxonList;
import dr.evomodel.tree.TreeModel;
import dr.xml.*;

import java.util.logging.Logger;

/**
 * @author Andrew Rambaut
 */
public class ExternalInternalBranchModelParser extends AbstractXMLObjectParser {

    public static final String EXTERNAL_INTERNAL_BRANCH_MODEL = "externalInternalBranchModel";
    public static final String EXTERNAL_BRANCHES = "externalBranches";

    public String getParserName() {
        return EXTERNAL_INTERNAL_BRANCH_MODEL;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        Logger.getLogger("dr.evomodel").info("Using external-internal branch model.");

        TreeModel tree = (TreeModel) xo.getChild(TreeModel.class);
        SubstitutionModel internalSubstitutionModel = (SubstitutionModel) xo.getChild(SubstitutionModel.class);
        SubstitutionModel externalSubstitutionModel = (SubstitutionModel) xo.getElementFirstChild(EXTERNAL_BRANCHES);

        return new ExternalInternalBranchModel(tree, externalSubstitutionModel, internalSubstitutionModel);
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return
                "This element provides a branch model which allows different substitution models" +
                        "on internal and external branches of the tree.";
    }

    public Class getReturnType() {
        return BranchSpecificBranchModel.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(TreeModel.class, "The tree"),
            new ElementRule(SubstitutionModel.class, "The substitution model for internal branches"),
            new ElementRule(EXTERNAL_BRANCHES,
                    new XMLSyntaxRule[]{
                            new ElementRule(SubstitutionModel.class, "The external branch substitution model"),
                    }, false)
    };

}
