/*
 * MATDataLikelihoodParser.java
 *
 * Copyright © 2026, the BEAST Development Team.
 * http://beast.community/about
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

package dr.evomodelxml.treedatalikelihood.MAT;

import dr.evolution.alignment.PatternList;
import dr.evomodel.branchmodel.BranchModel;
import dr.evomodel.branchmodel.HomogeneousBranchModel;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.siteratemodel.SiteRateModel;
import dr.evomodel.substmodel.SubstitutionModel;
import dr.evomodel.tree.TreeModel;
import dr.evomodel.treedatalikelihood.MAT.BranchHistoryModel;
import dr.evomodel.treedatalikelihood.MAT.ExplicitHitsHistoryModel;
import dr.evomodel.treedatalikelihood.MAT.IntegratedHitsHistoryModel;
import dr.evomodel.treedatalikelihood.MAT.MATDataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.MAT.SingleHitHistoryModel;
import dr.evomodel.treedatalikelihood.MAT.TimedHitsHistoryModel;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.inference.model.Parameter;
import dr.xml.*;

/**
 * XML parser for a TreeDataLikelihood using the mutation annotated tree delegate.
 */
public class MATDataLikelihoodParser extends AbstractXMLObjectParser {

    public static final String MAT_DATA_LIKELIHOOD = "matDataLikelihood";
    public static final String SITE_RATES = "siteRates";
    public static final String BRANCH_HISTORY = "branchHistory";

    @Override
    public String getParserName() {
        return MAT_DATA_LIKELIHOOD;
    }

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        PatternList patternList = (PatternList) xo.getChild(PatternList.class);
        TreeModel treeModel = (TreeModel) xo.getChild(TreeModel.class);

        BranchModel branchModel = (BranchModel) xo.getChild(BranchModel.class);
        if (branchModel == null) {
            SubstitutionModel substitutionModel = (SubstitutionModel) xo.getChild(SubstitutionModel.class);
            if (substitutionModel == null) {
                throw new XMLParseException("Either a branchModel or a substitution model must be supplied for " + xo.getId());
            }
            branchModel = new HomogeneousBranchModel(substitutionModel);
        }

        SiteRateModel siteRateModel = (SiteRateModel) xo.getChild(SiteRateModel.class);
        BranchRateModel branchRateModel = (BranchRateModel) xo.getChild(BranchRateModel.class);
        Parameter siteRates = xo.hasChildNamed(SITE_RATES) ? (Parameter) xo.getElementFirstChild(SITE_RATES) : null;

        final String historyName = xo.getAttribute(BRANCH_HISTORY, IntegratedHitsHistoryModel.NAME);
        if (historyName.equalsIgnoreCase(IntegratedHitsHistoryModel.NAME)) {
            final boolean rateCategories = siteRateModel != null && siteRateModel.getCategoryCount() > 1;
            if (siteRates != null || rateCategories) {
                throw new XMLParseException("Site rate heterogeneity (" +
                        (siteRates != null ? "a per-site rate parameter" :
                                "a site model with " + siteRateModel.getCategoryCount() + " rate categories") +
                        ") is not yet implemented for the " + IntegratedHitsHistoryModel.NAME + " branch history model" +
                        (xo.hasAttribute(BRANCH_HISTORY) ? "" : " (the default)") + " in " + xo.getId() +
                        ": use " + BRANCH_HISTORY + "=\"" + SingleHitHistoryModel.NAME + "\", \"" +
                        ExplicitHitsHistoryModel.NAME + "\" or \"" + TimedHitsHistoryModel.NAME +
                        "\" with a per-site rate parameter (" + SITE_RATES + ")");
            }
        }

        MATDataLikelihoodDelegate delegate;
        try {
            BranchHistoryModel historyModel = BranchHistoryModel.create(historyName,
                    patternList.getDataType().getStateCount());
            delegate = new MATDataLikelihoodDelegate(treeModel, patternList, branchModel,
                    siteRateModel, branchRateModel, siteRates, historyModel);
        } catch (IllegalArgumentException iae) {
            throw new XMLParseException(iae.getMessage() + " in " + xo.getId());
        }

        return new TreeDataLikelihood(delegate, treeModel, delegate.getBranchRateModel());
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return new XMLSyntaxRule[] {
                AttributeRule.newStringRule(BRANCH_HISTORY, true, "How each site's history along a branch is " +
                        "represented: " + String.join(", ", BranchHistoryModel.NAMES) + " (default " +
                        IntegratedHitsHistoryModel.NAME + ", which does not yet support site rate heterogeneity)"),
                new ElementRule(PatternList.class),
                new ElementRule(TreeModel.class),
                new XORRule(
                        new ElementRule(BranchModel.class),
                        new ElementRule(SubstitutionModel.class)
                ),
                new ElementRule(SiteRateModel.class, true),
                new ElementRule(BranchRateModel.class, true),
                new ElementRule(SITE_RATES, Parameter.class,
                        "Per-pattern relative rates (one per pattern), multiplying the site rate model rate", true)
        };
    }

    @Override
    public String getParserDescription() {
        return "Constructs a TreeDataLikelihood whose likelihood is that of a parsimonious mutation annotated tree.";
    }

    @Override
    @SuppressWarnings("rawtypes")
    public Class getReturnType() {
        return TreeDataLikelihood.class;
    }
}
