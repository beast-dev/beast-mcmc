/*
 * SequenceDistanceStatisticParser.java
 *
 * Copyright (c) 2002-2021 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.evomodelxml.treedatalikelihood;

import dr.evolution.alignment.PatternList;
import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeUtils;
import dr.evolution.util.Taxa;
import dr.evolution.util.TaxonList;
import dr.evomodel.substmodel.SubstitutionModel;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.evomodel.treedatalikelihood.discrete.SequenceDistanceStatistic;
import dr.evomodel.treelikelihood.AncestralStateBeagleTreeLikelihood;
import dr.inference.model.Statistic;
import dr.oldevomodel.sitemodel.SiteModel;
import dr.xml.*;

/**
 */
public class SequenceDistanceStatisticParser extends AbstractXMLObjectParser {

    public static String REPORT_DISTANCE = "reportDistance";
    public static String SEQUENCE_DISTANCE_STATISTIC = "sequenceDistanceStatistic";
    public static final String TAXA = "taxa";
    private static final String MRCA = "mrca";
    public static final String TREE_SEQUENCE_IS_ANCESTRAL = "treeSequenceIsAncestral";

    public String getParserName() { return SEQUENCE_DISTANCE_STATISTIC; }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

//        String name = xo.getAttribute(Statistic.NAME, xo.getId());

//        System.err.println("Got name = " + name);

        AncestralStateBeagleTreeLikelihood asrLike = (AncestralStateBeagleTreeLikelihood) xo.getChild(AncestralStateBeagleTreeLikelihood.class);

        SubstitutionModel subsModel = (SubstitutionModel) xo.getChild(SubstitutionModel.class);

//        if ( asrLike.getModelCount() != 1 && subsModel == null ) {
//            throw new RuntimeException("There are " + asrLike.getBranchRateModel().getModelName() + " substitution models. Partitioning not implemented, if using substitution model for tree data dataset must be unpartitioned.");
//        }

        PatternList patternList = (PatternList)xo.getChild(PatternList.class);

        // If true, sequence at given tree node is taken to be ancestral to user-supplied sequence
        // Otherwise, user-defined sequence is taken to be ancestral to sequence at tree node
        boolean treeSequenceIsAncestral = xo.getAttribute(TREE_SEQUENCE_IS_ANCESTRAL, false);

        // If true, distance between node and sequence are reported, if false the maximized likelihood
//        boolean reportDistance = xo.getAttribute(REPORT_DISTANCE, true);
        SequenceDistanceStatistic.DistanceType type = parseFromString(xo.getAttribute(REPORT_DISTANCE,
                SequenceDistanceStatistic.DistanceType.MAXIMIZED_DISTANCE.getName()));

        TaxonList taxa = null;
        if (xo.hasChildNamed(MRCA)) {
            taxa = (TaxonList) xo.getElementFirstChild(MRCA);
        }

        SequenceDistanceStatistic seqDistStatistic = null;
        try {
            seqDistStatistic = new SequenceDistanceStatistic(asrLike,subsModel,
                    patternList,treeSequenceIsAncestral, taxa, type);
        } catch (TreeUtils.MissingTaxonException e) {
            throw new XMLParseException("Unable to find taxon-set");
        }

        return seqDistStatistic;
    }

    public  SequenceDistanceStatistic.DistanceType parseFromString(String text) throws XMLParseException {
        for (SequenceDistanceStatistic.DistanceType type : SequenceDistanceStatistic.DistanceType.values()) {
            if (type.getName().toLowerCase().compareToIgnoreCase(text) == 0) {
                return type;
            }
        }
        throw new XMLParseException("Unknown type '" + text + "'");
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "A statistic .";
    }

    public Class getReturnType() { return SequenceDistanceStatistic.class; }

    public XMLSyntaxRule[] getSyntaxRules() { return rules; }

    private XMLSyntaxRule[] rules = new XMLSyntaxRule[] {
//            new StringAttributeRule(Statistic.NAME, "A name for this statistic primarily for the purposes of logging", true),
            new ElementRule(AncestralStateBeagleTreeLikelihood.class, false),
            new ElementRule(SubstitutionModel.class, true),
            new ElementRule(PatternList.class, false),
            AttributeRule.newBooleanRule(TREE_SEQUENCE_IS_ANCESTRAL, true),
            AttributeRule.newBooleanRule(REPORT_DISTANCE, true),
            new ElementRule(MRCA,
                    new XMLSyntaxRule[]{new ElementRule(Taxa.class)}, true),
    };

}
