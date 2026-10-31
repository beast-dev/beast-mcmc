/*
 * ALSTreeLikelihoodParser.java
 *
 * Copyright © 2002-2024 the BEAST Development Team
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
 *
 */

package dr.evomodelxml.treelikelihood;

import dr.evolution.alignment.PatternList;
import dr.evolution.util.Taxon;
import dr.evomodel.MSSD.ALSBeagleTreeLikelihood;
import dr.evomodel.MSSD.AbstractObservationProcess;
import dr.evomodel.MSSD.AnyTipObservationProcess;
import dr.evomodel.MSSD.SingleTipObservationProcess;
import dr.evomodel.branchmodel.BranchModel;
import dr.evomodel.branchmodel.HomogeneousBranchModel;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.siteratemodel.DiscretizedSiteRateModel;
import dr.evomodel.siteratemodel.SiteRateModel;
import dr.evomodel.substmodel.MutationDeathModel;
import dr.evomodel.substmodel.SubstitutionModel;
import dr.evomodel.tree.TreeModel;
import dr.evomodel.treelikelihood.PartialsRescalingScheme;
import dr.inference.model.Parameter;
import dr.xml.*;

import java.util.logging.Logger;

import static dr.evomodelxml.treelikelihood.BeagleTreeLikelihoodParser.USE_AMBIGUITIES;

/**
 *
 */
public class ALSTreeLikelihoodParser extends AbstractXMLObjectParser {
    public static final String LIKE_NAME = "alsTreeLikelihood";
    public static final String INTEGRATE_GAIN_RATE = "integrateGainRate";
    public static final String OBSERVATION_PROCESS = "observationProcess";
    public static final String OBSERVATION_TYPE = "type";
    public static final String OBSERVATION_TAXON = "taxon";
    public static final String ANY_TIP = "anyTip";
    public final static String IMMIGRATION_RATE = "immigrationRate";
    public static final String FORCE_RESCALING = "forceRescaling";
    public static final String STORE_PARTIALS = "storePartials";


    public String getParserName() {
        return LIKE_NAME;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        boolean useAmbiguities = false;
        boolean storePartials = true;
        if (xo.hasAttribute(USE_AMBIGUITIES)) {
            useAmbiguities = xo.getBooleanAttribute(USE_AMBIGUITIES);
        }
        if (xo.hasAttribute(STORE_PARTIALS)) {
            storePartials = xo.getBooleanAttribute(STORE_PARTIALS);
        }

        boolean integrateGainRate = xo.getBooleanAttribute(INTEGRATE_GAIN_RATE);

        //AbstractObservationProcess observationProcess = (AbstractObservationProcess) xo.getChild(AbstractObservationProcess.class);


        PatternList patternList = (PatternList) xo.getChild(PatternList.class);
        TreeModel treeModel = (TreeModel) xo.getChild(TreeModel.class);
        DiscretizedSiteRateModel siteModel = (DiscretizedSiteRateModel) xo.getChild(DiscretizedSiteRateModel.class);
        BranchModel branchModel = (BranchModel) xo.getChild(BranchModel.class);
        if (branchModel == null) {
            SubstitutionModel substitutionModel = (SubstitutionModel) xo.getChild(SubstitutionModel.class);
            if (substitutionModel == null) {
                substitutionModel = siteModel.getSubstitutionModel();
            }
            if (substitutionModel == null) {
                throw new XMLParseException("No substitution model available for TreeLikelihood: "+xo.getId());
            }
            branchModel = new HomogeneousBranchModel(substitutionModel);
        }
        BranchRateModel branchRateModel = (BranchRateModel) xo.getChild(BranchRateModel.class);
        Parameter mu = ((MutationDeathModel) siteModel.getSubstitutionModel()).getDeathParameter();
        Parameter lam;
        if (!integrateGainRate) {
            lam = (Parameter) xo.getElementFirstChild(IMMIGRATION_RATE);
        } else {
            lam = new Parameter.Default("gainRate", 1.0, 0.001, 1.999);
        }
        AbstractObservationProcess observationProcess = null;

        Logger.getLogger("dr.evolution").info("\n ---------------------------------\nCreating ALSTreeLikelihood model.");
        for (int i = 0; i < xo.getChildCount(); ++i) {
            Object cxo = xo.getChild(i);
            if (cxo instanceof XMLObject && ((XMLObject) cxo).getName().equals(OBSERVATION_PROCESS)) {
                if (((XMLObject) cxo).getStringAttribute(OBSERVATION_TYPE).equals("singleTip")) {
                    String taxonName = ((XMLObject) cxo).getStringAttribute(OBSERVATION_TAXON);
                    Taxon taxon = treeModel.getTaxon(treeModel.getTaxonIndex(taxonName));
                    observationProcess = new SingleTipObservationProcess(treeModel, patternList, siteModel,
                            branchRateModel, mu, lam, taxon);
                    Logger.getLogger("dr.evolution").info("All traits are assumed extant in " + taxonName);
                } else {  // "anyTip" observation process
                    observationProcess = new AnyTipObservationProcess(ANY_TIP, treeModel, patternList,
                            siteModel, branchRateModel, mu, lam);
                    Logger.getLogger("dr.evolution").info("Observed traits are assumed to be extant in at least one tip node.");
                }

                observationProcess.setIntegrateGainRate(integrateGainRate);
            }
        }
        Logger.getLogger("dr.evolution").info("\tIf you publish results using Acquisition-Loss-Mutation (ALS) Model likelihood, please reference Alekseyenko, Lee and Suchard (2008) Syst. Biol 57: 772-784.\n---------------------------------\n");

        PartialsRescalingScheme rescalingScheme = PartialsRescalingScheme.DEFAULT;
        if (xo.getAttribute(FORCE_RESCALING, false)) {
            rescalingScheme = PartialsRescalingScheme.ALWAYS;

        }
        return new ALSBeagleTreeLikelihood(observationProcess, patternList, treeModel, branchModel,
                siteModel, branchRateModel, null, useAmbiguities, rescalingScheme, true, null);
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "This element represents the likelihood of a patternlist on a tree given the site model.";
    }

    public Class getReturnType() {
        return ALSBeagleTreeLikelihood.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            AttributeRule.newBooleanRule(USE_AMBIGUITIES, true),
            AttributeRule.newBooleanRule(STORE_PARTIALS, true),
            AttributeRule.newBooleanRule(INTEGRATE_GAIN_RATE),
            AttributeRule.newBooleanRule(FORCE_RESCALING, true),
            new ElementRule(IMMIGRATION_RATE, new XMLSyntaxRule[]{new ElementRule(Parameter.class)}, true),
            new ElementRule(PatternList.class),
            new ElementRule(TreeModel.class),
            new ElementRule(DiscretizedSiteRateModel.class),
            new ElementRule(BranchModel.class, true),
            new ElementRule(BranchRateModel.class, true),
            new ElementRule(OBSERVATION_PROCESS,
                    new XMLSyntaxRule[]{AttributeRule.newStringRule(OBSERVATION_TYPE, false),
                            AttributeRule.newStringRule(OBSERVATION_TAXON, true)})
            //new ElementRule(AbstractObservationProcess.class)
    };

}
