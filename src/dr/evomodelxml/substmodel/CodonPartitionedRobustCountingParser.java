/*
 * CodonPartitionedRobustCountingParser.java
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

package dr.evomodelxml.substmodel;

import dr.evomodel.substmodel.CodonLabeling;
import dr.evomodel.substmodel.CodonPartitionedRobustCounting;
import dr.evomodel.substmodel.StratifiedTraitOutputFormat;
import dr.evomodel.treelikelihood.AncestralStateBeagleTreeLikelihood;
import dr.evolution.datatype.Codons;
import dr.evolution.datatype.GeneticCode;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.tree.TreeModel;
import dr.xml.*;

/**
 * @author Marc A. Suchard
 */
public class CodonPartitionedRobustCountingParser extends AbstractXMLObjectParser {

    public static final String PARSER_NAME = "codonPartitionedRobustCounting";
    public static final String FIRST = "firstPosition";
    public static final String SECOND = "secondPosition";
    public static final String THIRD = "thirdPosition";
    public static final String LABELING = "labeling";
    public static final String BRANCH_FORMAT = "branchFormat";
    public static final String LOG_FORMAT = "logFormat";
    public static final String USE_UNIFORMIZATION = "useUniformization";
    public static final String INCLUDE_EXTERNAL = "includeExternalBranches";
    public static final String INCLUDE_INTERNAL = "includeInternalBranches";
    public static final String DO_UNCONDITIONED_PER_BRANCH = "unconditionedPerBranch";
    public static final String SAVE_HISTORY = "saveCompleteHistory";
    public static final String AVERAGE_RATES = "averageRates";
    public static final String USE_NEW_NEUTRAL_MODEL = "useNewNeutralModel";
    public static final String PREFIX = "prefix";

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        AncestralStateBeagleTreeLikelihood[] partition = new AncestralStateBeagleTreeLikelihood[3];
        String[] labels = new String[]{FIRST, SECOND, THIRD};

        int patternCount = -1;

        BranchRateModel testBranchRateModel = null;
        for (int i = 0; i < 3; i++) {
            partition[i] = (AncestralStateBeagleTreeLikelihood)
                    xo.getChild(labels[i]).getChild(AncestralStateBeagleTreeLikelihood.class);

            if (i == 0) {
                patternCount = partition[i].getPatternCount();
            } else {
                if (partition[i].getPatternCount() != patternCount) {
                    throw new XMLParseException("Codon-partitioned robust counting requires all partitions to have the same length." +
                            " Make sure that partitions include all unique sites and do not strip gaps.");
                }
            }
            // Ensure that siteRateModel has one category
            if (partition[i].getSiteRateModel().getCategoryCount() > 1) {
                throw new XMLParseException("Robust counting currently only implemented for single category models");
            }

            // Ensure that branchRateModel is the same across all partitions
            if (testBranchRateModel == null) {
                testBranchRateModel = partition[i].getBranchRateModel();
            } else if (testBranchRateModel != partition[i].getBranchRateModel()) {
                throw new XMLParseException(
                        "Robust counting currently requires the same branch rate model for all partitions");
            }
        }

        TreeModel tree = (TreeModel) xo.getChild(TreeModel.class);

        Codons codons = Codons.UNIVERSAL;
        if (xo.hasAttribute(GeneticCode.GENETIC_CODE)) {
            String codeStr = xo.getStringAttribute(GeneticCode.GENETIC_CODE);
            codons = Codons.findByName(codeStr);
        }

        String labelingString = (String) xo.getAttribute(LABELING);
        CodonLabeling codonLabeling = CodonLabeling.parseFromString(labelingString);
        if (codonLabeling == null) {
            throw new XMLParseException("Unrecognized codon labeling '" + labelingString + "'");
        }

        String branchFormatString = xo.getAttribute(BRANCH_FORMAT,
                StratifiedTraitOutputFormat.SUM_OVER_SITES.getText());
        StratifiedTraitOutputFormat branchFormat = StratifiedTraitOutputFormat.parseFromString(branchFormatString);
        if (branchFormat == null) {
            throw new XMLParseException("Unrecognized branch output format '" + branchFormat + "'");
        }

        String logFormatString = xo.getAttribute(LOG_FORMAT,
                StratifiedTraitOutputFormat.SUM_OVER_SITES.getText());
        StratifiedTraitOutputFormat logFormat = StratifiedTraitOutputFormat.parseFromString(logFormatString);
        if (logFormat == null) {
            throw new XMLParseException("Unrecognized log output format '" + branchFormat + "'");
        }

        boolean useUniformization = xo.getAttribute(USE_UNIFORMIZATION, false);
        boolean includeExternalBranches = xo.getAttribute(INCLUDE_EXTERNAL, true);
        boolean includeInternalBranches = xo.getAttribute(INCLUDE_INTERNAL, true);
        boolean doUnconditionedPerBranch = xo.getAttribute(DO_UNCONDITIONED_PER_BRANCH, false);
        boolean averageRates = xo.getAttribute(AVERAGE_RATES, true);
        boolean saveCompleteHistory = xo.getAttribute(SAVE_HISTORY, false);
        boolean useNewNeutralModel = xo.getAttribute(USE_NEW_NEUTRAL_MODEL, false);

        String prefix = xo.hasAttribute(PREFIX) ? xo.getStringAttribute(PREFIX) : null;

        return new CodonPartitionedRobustCounting(
                xo.getId(),
                tree,
                partition,
                codons,
                codonLabeling,
                useUniformization,
                includeExternalBranches,
                includeInternalBranches,
                doUnconditionedPerBranch,
                saveCompleteHistory,
                averageRates,
                useNewNeutralModel,
                branchFormat,
                logFormat,
                prefix);
    }

    private static final XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
            new ElementRule(FIRST,
                    new XMLSyntaxRule[]{
                            new ElementRule(AncestralStateBeagleTreeLikelihood.class),
                    }),
            new ElementRule(SECOND,
                    new XMLSyntaxRule[]{
                            new ElementRule(AncestralStateBeagleTreeLikelihood.class),
                    }),
            new ElementRule(THIRD,
                    new XMLSyntaxRule[]{
                            new ElementRule(AncestralStateBeagleTreeLikelihood.class),
                    }),
            new ElementRule(TreeModel.class),
            new StringAttributeRule(GeneticCode.GENETIC_CODE,
                    "The genetic code to use",
                    GeneticCode.GENETIC_CODE_NAMES, true),
            AttributeRule.newBooleanRule(USE_UNIFORMIZATION, true),
            AttributeRule.newBooleanRule(INCLUDE_EXTERNAL, true),
            AttributeRule.newBooleanRule(INCLUDE_INTERNAL, true),
            AttributeRule.newBooleanRule(DO_UNCONDITIONED_PER_BRANCH, true),
            AttributeRule.newBooleanRule(AVERAGE_RATES, true),
            AttributeRule.newStringRule(LABELING),
            AttributeRule.newBooleanRule(SAVE_HISTORY, true),
            AttributeRule.newBooleanRule(USE_NEW_NEUTRAL_MODEL, true),
    };

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    public String getParserDescription() {
        return "A parser to specify robust counting procedures on codon partitioned models";
    }

    public Class getReturnType() {
        return CodonPartitionedRobustCounting.class;
    }

    public String getParserName() {
        return PARSER_NAME;
    }
}
