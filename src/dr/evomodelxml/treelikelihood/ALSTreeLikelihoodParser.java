/*
 * ALSTreeLikelihoodParser.java
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

package dr.evomodelxml.treelikelihood;

import dr.evolution.alignment.PatternList;
import dr.evolution.util.Taxon;
import dr.oldevomodel.MSSD.ALSTreeLikelihood;
import dr.oldevomodel.MSSD.AbstractObservationProcess;
import dr.oldevomodel.MSSD.AnyTipObservationProcess;
import dr.oldevomodel.MSSD.SingleTipObservationProcess;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.oldevomodel.sitemodel.SiteModel;
import dr.oldevomodel.substmodel.MutationDeathModel;
import dr.evomodel.tree.TreeModel;
import dr.oldevomodelxml.treelikelihood.TreeLikelihoodParser;
import dr.inference.model.Parameter;
import dr.xml.*;

import java.util.logging.Logger;

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
    public static final String FORCE_RESCALING = TreeLikelihoodParser.FORCE_RESCALING;

    public String getParserName() {
        return LIKE_NAME;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        boolean useAmbiguities = false;
        boolean storePartials = true;
        if (xo.hasAttribute(TreeLikelihoodParser.USE_AMBIGUITIES)) {
            useAmbiguities = xo.getBooleanAttribute(TreeLikelihoodParser.USE_AMBIGUITIES);
        }
        if (xo.hasAttribute(TreeLikelihoodParser.STORE_PARTIALS)) {
            storePartials = xo.getBooleanAttribute(TreeLikelihoodParser.STORE_PARTIALS);
        }

        boolean integrateGainRate = xo.getBooleanAttribute(INTEGRATE_GAIN_RATE);

        //AbstractObservationProcess observationProcess = (AbstractObservationProcess) xo.getChild(AbstractObservationProcess.class);


        PatternList patternList = (PatternList) xo.getChild(PatternList.class);
        TreeModel treeModel = (TreeModel) xo.getChild(TreeModel.class);
        SiteModel siteModel = (SiteModel) xo.getChild(SiteModel.class);
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

        boolean forceRescaling = xo.getAttribute(FORCE_RESCALING, false);

//        forceRescaling = true;

        return new ALSTreeLikelihood(observationProcess, patternList, treeModel, siteModel, branchRateModel,
                useAmbiguities, storePartials, forceRescaling);
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "This element represents the likelihood of a patternlist on a tree given the site model.";
    }

    public Class getReturnType() {
        return ALSTreeLikelihood.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            AttributeRule.newBooleanRule(TreeLikelihoodParser.USE_AMBIGUITIES, true),
            AttributeRule.newBooleanRule(TreeLikelihoodParser.STORE_PARTIALS, true),
            AttributeRule.newBooleanRule(INTEGRATE_GAIN_RATE),
            AttributeRule.newBooleanRule(FORCE_RESCALING, true),
            new ElementRule(IMMIGRATION_RATE, new XMLSyntaxRule[]{new ElementRule(Parameter.class)}, true),
            new ElementRule(PatternList.class),
            new ElementRule(TreeModel.class),
            new ElementRule(SiteModel.class),
            new ElementRule(BranchRateModel.class, true),
            new ElementRule(OBSERVATION_PROCESS,
                    new XMLSyntaxRule[]{AttributeRule.newStringRule(OBSERVATION_TYPE, false),
                            AttributeRule.newStringRule(OBSERVATION_TAXON, true)})
            //new ElementRule(AbstractObservationProcess.class)
    };

}
