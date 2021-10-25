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
import dr.evolution.datatype.Nucleotides;
import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeUtils;
import dr.evolution.util.Taxa;
import dr.evolution.util.TaxonList;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.branchratemodel.StrictClockBranchRates;
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
    private static final String MRCA = "mrca";
    public static final String TAXA = "taxa";

    public String getParserName() { return SEQUENCE_DISTANCE_STATISTIC; }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        AncestralStateBeagleTreeLikelihood asrLike = (AncestralStateBeagleTreeLikelihood) xo.getChild(AncestralStateBeagleTreeLikelihood.class);

        SubstitutionModel subsModel = (SubstitutionModel) xo.getChild(SubstitutionModel.class);

        BranchRateModel branchRates = (BranchRateModel)xo.getChild(BranchRateModel.class);

        if ( !(branchRates instanceof  StrictClockBranchRates)) {
            throw new XMLParseException("Clock models other than StrictClockBranchRates not currently supported.");
        }

        PatternList patternList = (PatternList)xo.getChild(PatternList.class);

        if ( patternList.areUnique() ) {
            throw new XMLParseException("Sequences being compared to tree nodes cannot be compressed (unique) patterns.");
        }

        // If true, distance between node and sequence are reported, if false the maximized likelihood
        SequenceDistanceStatistic.DistanceType type = parseFromString(xo.getAttribute(REPORT_DISTANCE,
                SequenceDistanceStatistic.DistanceType.MAXIMIZED_DISTANCE.getName()));

        TaxonList mrcaTaxa = null;
        if (xo.hasChildNamed(MRCA)) {
            mrcaTaxa = (TaxonList) xo.getElementFirstChild(MRCA);
        }

        SequenceDistanceStatistic seqDistStatistic = null;
        try {
            seqDistStatistic = new SequenceDistanceStatistic(asrLike,subsModel,branchRates,
                    patternList, mrcaTaxa, type);
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
            new ElementRule(AncestralStateBeagleTreeLikelihood.class, false),
            new ElementRule(SubstitutionModel.class, true),
            new ElementRule(PatternList.class, true),
            new ElementRule(BranchRateModel.class, false),
            AttributeRule.newBooleanRule(REPORT_DISTANCE, true),
            new ElementRule(MRCA,
                    new XMLSyntaxRule[]{new ElementRule(Taxa.class)}, true),
    };

}
