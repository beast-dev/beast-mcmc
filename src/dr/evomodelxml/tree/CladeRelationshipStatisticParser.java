/*
 * CladeRelationStatisticParser.java
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

package dr.evomodelxml.tree;

import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeUtils;
import dr.evolution.util.Taxa;
import dr.evolution.util.Taxon;
import dr.evolution.util.TaxonList;
import dr.evomodel.speciation.BirthDeathGernhard08Model;
import dr.evomodel.speciation.NewBirthDeathSerialSamplingModel;
import dr.evomodel.speciation.SpeciationLikelihoodGradient;
import dr.evomodel.tree.AncestorOnStemStatistic;
import dr.evomodel.tree.AncestralTraitTreeModel;
import dr.evomodel.tree.CladeRelationshipStatistic;
import dr.inference.model.Statistic;
import dr.xml.*;

import static dr.evomodelxml.tree.MonophylyStatisticParser.parseTaxonListOrTaxa;

/**
 */
public class CladeRelationshipStatisticParser extends AbstractXMLObjectParser {

    public static final String ANCESTOR_STATISTIC = "cladeRelationshipStatistic";
    private static final String MRCA_A = "taxaA";
    private static final String MRCA_B = "taxaB";
    private static final String RELATIONSHIP_TYPE = "relationshipType";

    public String getParserName() {
        return ANCESTOR_STATISTIC;
    }


    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        String name = xo.getAttribute(Statistic.NAME, xo.getId());

        Tree tree = (Tree) xo.getChild(Tree.class);

        TaxonList taxaA = null;
        if (xo.hasChildNamed(MRCA_A)) {
            taxaA = parseTaxonListOrTaxa(xo.getChild(MRCA_A));
        }

        TaxonList taxaB = null;
        if (xo.hasChildNamed(MRCA_B)) {
            taxaB = parseTaxonListOrTaxa(xo.getChild(MRCA_B));
        }

        String relationshipType = (String) xo.getAttribute(RELATIONSHIP_TYPE);

        CladeRelationshipStatistic.RelationshipType type = CladeRelationshipStatistic.factory(relationshipType);


        try {
            return new CladeRelationshipStatistic(name, tree, taxaA, taxaB, type);
        } catch (TreeUtils.MissingTaxonException mte) {
            throw new XMLParseException("Taxon, " + mte + ", in " + getParserName() + "was not found in the tree.");
        }
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "A statistic that returns true if the clades defined by taxon sets A and B have the specified relationship.";
    }

    public Class getReturnType() {
        return CladeRelationshipStatistic.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new StringAttributeRule(Statistic.NAME, "A name for this statistic for the purpose of logging", true),
            new ElementRule(Tree.class),
            new ElementRule(MRCA_A, new XMLSyntaxRule[]{
                    new XORRule(
                            new ElementRule(Taxon.class, 1, Integer.MAX_VALUE),
                            new ElementRule(Taxa.class)
                    )
            }),
            new ElementRule(MRCA_B, new XMLSyntaxRule[]{
                    new XORRule(
                            new ElementRule(Taxon.class, 1, Integer.MAX_VALUE),
                            new ElementRule(Taxa.class)
                    )
            })
    };

}
