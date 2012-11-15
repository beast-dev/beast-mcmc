/*
 * TreeLikelihoodParser.java
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

//import dr.app.beagle.evomodel.treelikelihood.RestrictedPartialsSequenceLikelihood;

import dr.app.beagle.evomodel.sitemodel.BranchSubstitutionModel;
import dr.app.beagle.evomodel.sitemodel.GammaSiteRateModel;
import dr.app.beagle.evomodel.sitemodel.HomogenousBranchSubstitutionModel;
import dr.app.beagle.evomodel.substmodel.FrequencyModel;
import dr.app.beagle.evomodel.treelikelihood.OldBeagleTreeLikelihood;
import dr.app.beagle.evomodel.treelikelihood.PartialsRescalingScheme;
import dr.evolution.alignment.PatternList;
import dr.evolution.alignment.Patterns;
import dr.evolution.alignment.SitePatterns;
import dr.evolution.tree.Tree;
import dr.evolution.util.TaxonList;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.tree.TreeModel;
import dr.evomodel.treelikelihood.TipStatesModel;
import dr.inference.model.CompoundLikelihood;
import dr.inference.model.Likelihood;
import dr.inference.model.Parameter;
import dr.xml.*;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.Set;

/**
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @author Marc Suchard
 * @version $Id$
 */

@Deprecated // Switching to BranchModel
public class OldTreeLikelihoodParser extends AbstractXMLObjectParser {

    public static final String BEAGLE_INSTANCE_COUNT = "beagle.instance.count";

    public static final String TREE_LIKELIHOOD = "oldTreeLikelihood";
    public static final String USE_AMBIGUITIES = "useAmbiguities";
    public static final String INSTANCE_COUNT = "instanceCount";
    //    public static final String DEVICE_NUMBER = "deviceNumber";
//    public static final String PREFER_SINGLE_PRECISION = "preferSinglePrecision";
    public static final String SCALING_SCHEME = "scalingScheme";
    public static final String PARTIALS_RESTRICTION = "partialsRestriction";

    public String getParserName() {
        return TREE_LIKELIHOOD;
    }

    protected OldBeagleTreeLikelihood createTreeLikelihood(PatternList patternList, TreeModel treeModel,
                                                        BranchSubstitutionModel branchSubstitutionModel, GammaSiteRateModel siteRateModel,
                                                        BranchRateModel branchRateModel,
                                                        TipStatesModel tipStatesModel,
                                                        boolean useAmbiguities, PartialsRescalingScheme scalingScheme,
                                                        Map<Set<String>, Parameter> partialsRestrictions,
                                                        XMLObject xo) throws XMLParseException {
        return new OldBeagleTreeLikelihood(
                patternList,
                treeModel,
                branchSubstitutionModel,
                siteRateModel,
                branchRateModel,
                tipStatesModel,
                useAmbiguities,
                scalingScheme,
                partialsRestrictions
        );
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        boolean useAmbiguities = xo.getAttribute(USE_AMBIGUITIES, false);
        int instanceCount = xo.getAttribute(INSTANCE_COUNT, 1);
        if (instanceCount < 1) {
            instanceCount = 1;
        }

        String ic = System.getProperty(BEAGLE_INSTANCE_COUNT);
        if (ic != null && ic.length() > 0) {
            instanceCount = Integer.parseInt(ic);
        }

        PatternList patternList = (PatternList) xo.getChild(PatternList.class);
        TreeModel treeModel = (TreeModel) xo.getChild(TreeModel.class);
        GammaSiteRateModel siteRateModel = (GammaSiteRateModel) xo.getChild(GammaSiteRateModel.class);

        FrequencyModel rootFreqModel = (FrequencyModel) xo.getChild(FrequencyModel.class);

        BranchSubstitutionModel branchSubstitutionModel = (BranchSubstitutionModel) xo.getChild(BranchSubstitutionModel.class);
        if (branchSubstitutionModel == null) {
            branchSubstitutionModel = new HomogenousBranchSubstitutionModel(
                    siteRateModel.getSubstitutionModel(),
                    (rootFreqModel != null ? rootFreqModel :
                            siteRateModel.getSubstitutionModel().getFrequencyModel())
            );
        }

        BranchRateModel branchRateModel = (BranchRateModel) xo.getChild(BranchRateModel.class);

        TipStatesModel tipStatesModel = (TipStatesModel) xo.getChild(TipStatesModel.class);
//        if (xo.getChild(TipStatesModel.class) != null) {
//            throw new XMLParseException("Sequence Error Models are not supported under BEAGLE yet. Please use Native BEAST Likelihood.");
//        }

        PartialsRescalingScheme scalingScheme = PartialsRescalingScheme.DEFAULT;
        if (xo.hasAttribute(SCALING_SCHEME)) {
            scalingScheme = PartialsRescalingScheme.parseFromString(xo.getStringAttribute(SCALING_SCHEME));
            if (scalingScheme == null)
                throw new XMLParseException("Unknown scaling scheme '"+xo.getStringAttribute(SCALING_SCHEME)+"' in "+
                        "OldBeagleTreeLikelihood object '"+xo.getId());

        }

        Map<Set<String>, Parameter> partialsRestrictions = null;

        if (xo.hasChildNamed(PARTIALS_RESTRICTION)) {
            XMLObject cxo = xo.getChild(PARTIALS_RESTRICTION);
            TaxonList taxonList = (TaxonList) cxo.getChild(TaxonList.class);
//            Parameter parameter = (Parameter) cxo.getChild(Parameter.class);
            try {
                Tree.Utils.getLeavesForTaxa(treeModel, taxonList);
            } catch (Tree.MissingTaxonException e) {
                throw new XMLParseException("Unable to parse taxon list: " + e.getMessage());
            }
            throw new XMLParseException("Restricting internal nodes is not yet implemented.  Contact Marc");

        }


        if (instanceCount == 1 || patternList.getPatternCount() < instanceCount) {
            return createTreeLikelihood(
                    patternList,
                    treeModel,
                    branchSubstitutionModel,
                    siteRateModel,
                    branchRateModel,
                    tipStatesModel,
                    useAmbiguities,
                    scalingScheme,
                    partialsRestrictions,
                    xo
            );
        }

        // using multiple instances of BEAGLE...

        if (!(patternList instanceof SitePatterns)) {
            throw new XMLParseException("BEAGLE_INSTANCES option cannot be used with BEAUti-selected codon partitioning.");
        }

        if (tipStatesModel != null) {
            throw new XMLParseException("BEAGLE_INSTANCES option cannot be used with a TipStateModel (i.e., a sequence error model).");
        }

        List<Likelihood> likelihoods = new ArrayList<Likelihood>();
        for (int i = 0; i < instanceCount; i++) {

            Patterns subPatterns = new Patterns((SitePatterns)patternList, 0, 0, 1, i, instanceCount);

            OldBeagleTreeLikelihood treeLikelihood = createTreeLikelihood(
                    subPatterns,
                    treeModel,
                    branchSubstitutionModel,
                    siteRateModel,
                    branchRateModel,
                    null,
                    useAmbiguities,
                    scalingScheme,
                    partialsRestrictions,
                    xo);
            treeLikelihood.setId(xo.getId() + "_" + instanceCount);
            likelihoods.add(treeLikelihood);
        }

        return new CompoundLikelihood(likelihoods);
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "This element represents the likelihood of a patternlist on a tree given the site model.";
    }

    public Class getReturnType() {
        return Likelihood.class;
    }

    public static final XMLSyntaxRule[] rules = {
            AttributeRule.newBooleanRule(USE_AMBIGUITIES, true),
            new ElementRule(PatternList.class),
            new ElementRule(TreeModel.class),
            new ElementRule(GammaSiteRateModel.class),
            new ElementRule(BranchSubstitutionModel.class, true),
            new ElementRule(BranchRateModel.class, true),
            new ElementRule(TipStatesModel.class, true),
            AttributeRule.newStringRule(SCALING_SCHEME,true),
            new ElementRule(PARTIALS_RESTRICTION, new XMLSyntaxRule[] {
                    new ElementRule(TaxonList.class),
                    new ElementRule(Parameter.class),
            }, true),
            new ElementRule(FrequencyModel.class, true),
            new ElementRule(TipStatesModel.class, true)
    };

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }
}
