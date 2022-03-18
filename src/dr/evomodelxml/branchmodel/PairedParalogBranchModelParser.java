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

import dr.evolution.util.Taxa;
import dr.evomodel.branchmodel.PairedParalogBranchModel;
import dr.evomodel.substmodel.SubstitutionModel;
import dr.evomodel.substmodel.geneconversion.PairedParalogGeneConversionSubstitutionModel;
import dr.evomodel.tree.DuplicationTreeModel;
import dr.xml.*;

/**
 * @author Xiang Ji
 * @author Jeff Thorne
 * @author Marc A. Suchard
 */
public class PairedParalogBranchModelParser extends AbstractXMLObjectParser {

    private final String NAME = "pairedParalogBranchModel";
    private final String BASE = "base";
    private final String GENE_CONVERSION = "geneConversion";

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        SubstitutionModel baseModel = (SubstitutionModel) xo.getChild(BASE).getChild(SubstitutionModel.class);
        SubstitutionModel geneConversionModel = (SubstitutionModel) xo.getChild(GENE_CONVERSION).getChild(SubstitutionModel.class);

        DuplicationTreeModel treeModel = (DuplicationTreeModel) xo.getChild(DuplicationTreeModel.class);
        Taxa postDuplicationTaxa = (Taxa) xo.getChild(Taxa.class);

        return new PairedParalogBranchModel(NAME, baseModel, geneConversionModel, treeModel, postDuplicationTaxa);
    }


    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(BASE, PairedParalogGeneConversionSubstitutionModel.class, "Base substitution model", false),
            new ElementRule(GENE_CONVERSION, PairedParalogGeneConversionSubstitutionModel.class, "Gene conversion substitution model", false),
            new ElementRule(Taxa.class, "All post duplication taxa (including duplication taxon)."),
            new ElementRule(DuplicationTreeModel.class, "An ancestralTraitTreeModel that tracks duplication node"),
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
