/*
 * BalancedBeagleTreeLikelihoodParser.java
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

package dr.evomodelxml.treelikelihood;

import dr.evolution.tree.TreeUtils;
import dr.evomodel.branchmodel.BranchModel;
import dr.evomodel.branchmodel.HomogeneousBranchModel;
import dr.evomodel.siteratemodel.GammaSiteRateModel;
import dr.evomodel.substmodel.FrequencyModel;
import dr.evomodel.substmodel.SubstitutionModel;
import dr.evomodel.treelikelihood.AbstractTreeLikelihood;
import dr.evomodel.treelikelihood.BeagleTreeLikelihood;
import dr.evomodel.treelikelihood.PartialsRescalingScheme;
import dr.evolution.alignment.PatternList;
import dr.evolution.alignment.Patterns;
import dr.evolution.alignment.SitePatterns;
import dr.evolution.util.TaxonList;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.tree.TreeModel;
import dr.evomodel.tipstatesmodel.TipStatesModel;
import dr.inference.model.CompoundLikelihood;
import dr.inference.model.Likelihood;
import dr.inference.model.Parameter;
import dr.xml.*;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.Set;

/**
 * @author Guy Baele
 */
public class BalancedBeagleTreeLikelihoodParser extends AbstractXMLObjectParser {

    //public static final String BEAGLE_INSTANCE_COUNT = "beagle.instance.count";

    public static final String TREE_LIKELIHOOD = "balancedTreeLikelihood";
    public static final String INSTANCE_COUNT = "instanceCount";
    public static final String PARTIALS_RESTRICTION = "partialsRestriction";
    
    public final int TEST_RUNS = 100;
    public final double TEST_CUTOFF = 1.30;

    public String getParserName() {
        return TREE_LIKELIHOOD;
    }

    protected BeagleTreeLikelihood createTreeLikelihood(PatternList patternList, TreeModel treeModel,
                                                        BranchModel branchModel,
                                                        GammaSiteRateModel siteRateModel,
                                                        BranchRateModel branchRateModel,
                                                        TipStatesModel tipStatesModel,
                                                        boolean useAmbiguities, PartialsRescalingScheme scalingScheme,
                                                        boolean delayScaling,
                                                        Map<Set<String>, Parameter> partialsRestrictions,
                                                        XMLObject xo) throws XMLParseException {
        return new BeagleTreeLikelihood(
                patternList,
                treeModel,
                branchModel,
                siteRateModel,
                branchRateModel,
                tipStatesModel,
                useAmbiguities,
                scalingScheme,
                delayScaling,
                partialsRestrictions
        );
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        boolean useAmbiguities = xo.getAttribute(BeagleTreeLikelihoodParser.USE_AMBIGUITIES, false);
        /*int instanceCount = xo.getAttribute(INSTANCE_COUNT, 1);
        if (instanceCount < 1) {
            instanceCount = 1;
        }

        String ic = System.getProperty(BEAGLE_INSTANCE_COUNT);
        if (ic != null && ic.length() > 0) {
            instanceCount = Integer.parseInt(ic);
        }*/

        PatternList patternList = (PatternList) xo.getChild(PatternList.class);
        TreeModel treeModel = (TreeModel) xo.getChild(TreeModel.class);
        GammaSiteRateModel siteRateModel = (GammaSiteRateModel) xo.getChild(GammaSiteRateModel.class);

        FrequencyModel rootFreqModel = (FrequencyModel) xo.getChild(FrequencyModel.class);

        BranchModel branchModel = (BranchModel) xo.getChild(BranchModel.class);
        if (branchModel == null) {
            SubstitutionModel substitutionModel = (SubstitutionModel) xo.getChild(SubstitutionModel.class);
            if (substitutionModel == null) {
                substitutionModel = siteRateModel.getSubstitutionModel();
            }
            if (substitutionModel == null) {
                throw new XMLParseException("No substitution model available for TreeLikelihood: "+xo.getId());
            }
            branchModel = new HomogeneousBranchModel(substitutionModel, rootFreqModel);
        }

        BranchRateModel branchRateModel = (BranchRateModel) xo.getChild(BranchRateModel.class);

        TipStatesModel tipStatesModel = (TipStatesModel) xo.getChild(TipStatesModel.class);
//        if (xo.getChild(TipStatesModel.class) != null) {
//            throw new XMLParseException("Sequence Error Models are not supported under BEAGLE yet. Please use Native BEAST Likelihood.");
//        }

        PartialsRescalingScheme scalingScheme = PartialsRescalingScheme.DEFAULT;
        if (xo.hasAttribute(BeagleTreeLikelihoodParser.SCALING_SCHEME)) {
//            scalingScheme = PartialsRescalingScheme.parseFromString(xo.getStringAttribute(BeagleTreeLikelihoodParser.SCALING_SCHEME));
            if (scalingScheme == null)
                throw new XMLParseException("Unknown scaling scheme '"+xo.getStringAttribute(BeagleTreeLikelihoodParser.SCALING_SCHEME)+"' in "+
                        "OldBeagleTreeLikelihood object '"+xo.getId());

        }

        boolean delayScaling = true;
        Map<Set<String>, Parameter> partialsRestrictions = null;

        if (xo.hasChildNamed(PARTIALS_RESTRICTION)) {
            XMLObject cxo = xo.getChild(PARTIALS_RESTRICTION);
            TaxonList taxonList = (TaxonList) cxo.getChild(TaxonList.class);
//            Parameter parameter = (Parameter) cxo.getChild(Parameter.class);
            try {
                TreeUtils.getLeavesForTaxa(treeModel, taxonList);
            } catch (TreeUtils.MissingTaxonException e) {
                throw new XMLParseException("Unable to parse taxon list: " + e.getMessage());
            }
            throw new XMLParseException("Restricting internal nodes is not yet implemented.  Contact Marc");

        }

        /*if (instanceCount == 1 || patternList.getPatternCount() < instanceCount) {
            return createTreeLikelihood(
                    patternList,
                    treeModel,
                    branchModel,
                    siteRateModel,
                    branchRateModel,
                    tipStatesModel,
                    useAmbiguities,
                    scalingScheme,
                    partialsRestrictions,
                    xo
            );
        }*/
        
        //first run a test for instanceCount == 1
        System.err.println("\nTesting instanceCount == 1");
        Likelihood baseLikelihood = createTreeLikelihood(
                patternList,
                treeModel,
                branchModel,
                siteRateModel,
                branchRateModel,
                tipStatesModel,
                useAmbiguities,
                scalingScheme,
                delayScaling,
                partialsRestrictions,
                xo
        );
        double start = System.nanoTime();
        for (int i = 0; i < TEST_RUNS; i++) {
        	baseLikelihood.makeDirty();
        	baseLikelihood.getLogLikelihood();
        }
        double end = System.nanoTime();
        double baseResult = end - start;
        System.err.println("Evaluation took: " + baseResult);

        // using multiple instances of BEAGLE...

        if (!(patternList instanceof SitePatterns)) {
            throw new XMLParseException("BEAGLE_INSTANCES option cannot be used with BEAUti-selected codon partitioning.");
        }

        if (tipStatesModel != null) {
            throw new XMLParseException("BEAGLE_INSTANCES option cannot be used with a TipStateModel (i.e., a sequence error model).");
        }

        //List<Likelihood> likelihoods = new ArrayList<Likelihood>();
        List<Likelihood> likelihoods = null;
        CompoundLikelihood compound = null;
        
        int instanceCount = 2;
        boolean optimal = false;
        
        while (optimal == false) {
        	
        	System.err.println("\nCreating instanceCount == " + instanceCount);
        	likelihoods = new ArrayList<Likelihood>();

        	for (int i = 0; i < instanceCount; i++) {
        		Patterns subPatterns = new Patterns((SitePatterns)patternList, 0, 0, 1, i, instanceCount);
        		
        		AbstractTreeLikelihood treeLikelihood = createTreeLikelihood(
                        subPatterns,
                        treeModel,
                        branchModel,
                        siteRateModel,
                        branchRateModel,
                        null,
                        useAmbiguities,
                        scalingScheme,
                        delayScaling,
                        partialsRestrictions,
                        xo);
                treeLikelihood.setId(xo.getId() + "_" + instanceCount);
                likelihoods.add(treeLikelihood);
        	}
        	
        	//construct compoundLikelihood
        	compound = new CompoundLikelihood(instanceCount, likelihoods);
        	
        	//test timings 
        	System.err.println("\nTesting instanceCount == " + instanceCount);
        	start = System.nanoTime();
            for (int i = 0; i < TEST_RUNS; i++) {
            	compound.makeDirty();
            	compound.getLogLikelihood();
            }
            end = System.nanoTime();
            double newResult = end - start;
            System.err.println("Evaluation took: " + newResult);
            
            if (baseResult/newResult > TEST_CUTOFF) {
            	
            	instanceCount++;
            	baseResult = newResult;
            	
            } else {
        	
            	optimal = true;
            	instanceCount--;
            	
            	System.err.println("\nCreating final BeagleTreeLikelihood with instanceCount: " + instanceCount);
            	
            	likelihoods = new ArrayList<Likelihood>();

            	for (int i = 0; i < instanceCount; i++) {
            		Patterns subPatterns = new Patterns((SitePatterns)patternList, 0, 0, 1, i, instanceCount);
            		
            		AbstractTreeLikelihood treeLikelihood = createTreeLikelihood(
                            subPatterns,
                            treeModel,
                            branchModel,
                            siteRateModel,
                            branchRateModel,
                            null,
                            useAmbiguities,
                            scalingScheme,
                            delayScaling,
                            partialsRestrictions,
                            xo);
                    treeLikelihood.setId(xo.getId() + "_" + instanceCount);
                    likelihoods.add(treeLikelihood);
            	}
            	
            	//construct compoundLikelihood
            	compound = new CompoundLikelihood(instanceCount, likelihoods);
        	
            }
        	
        }
        
        return compound;
        
        /*for (int i = 0; i < instanceCount; i++) {

            Patterns subPatterns = new Patterns((SitePatterns)patternList, 0, 0, 1, i, instanceCount);

            AbstractTreeLikelihood treeLikelihood = createTreeLikelihood(
                    subPatterns,
                    treeModel,
                    branchModel,
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

        return new CompoundLikelihood(likelihoods);*/
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "This element represents the likelihood of a patternlist on a tree given the site model, with an automated detection of instanceCount.";
    }

    public Class getReturnType() {
        return Likelihood.class;
    }

    public static final XMLSyntaxRule[] rules = {
            AttributeRule.newBooleanRule(BeagleTreeLikelihoodParser.USE_AMBIGUITIES, true),
            new ElementRule(PatternList.class),
            new ElementRule(TreeModel.class),
            new ElementRule(GammaSiteRateModel.class),
            new ElementRule(BranchModel.class, true),
            new ElementRule(SubstitutionModel.class, true),
            new ElementRule(BranchRateModel.class, true),
            new ElementRule(TipStatesModel.class, true),
            AttributeRule.newStringRule(BeagleTreeLikelihoodParser.SCALING_SCHEME,true),
            AttributeRule.newBooleanRule(BeagleTreeLikelihoodParser.DELAY_SCALING,true),
            new ElementRule(PARTIALS_RESTRICTION, new XMLSyntaxRule[] {
                    new ElementRule(TaxonList.class),
                    new ElementRule(Parameter.class),
            }, true),
            new ElementRule(TipStatesModel.class, true),
            new ElementRule(FrequencyModel.class, true),
    };

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }
}
