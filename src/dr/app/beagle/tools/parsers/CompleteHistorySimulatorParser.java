/*
 * CompleteHistorySimulatorParser.java
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

package dr.app.beagle.tools.parsers;

import dr.evomodelxml.treelikelihood.MarkovJumpsTreeLikelihoodParser;
import dr.evomodel.siteratemodel.GammaSiteRateModel;
import dr.evomodel.substmodel.CodonLabeling;
import dr.app.beagle.tools.CompleteHistorySimulator;
import dr.evolution.alignment.Alignment;
import dr.evolution.datatype.Codons;
import dr.evolution.datatype.DataType;
import dr.evolution.tree.Tree;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.branchratemodel.DefaultBranchRateModel;
import dr.inference.markovjumps.MarkovJumpsType;
import dr.inference.model.Parameter;
import dr.xml.*;

/**
 * @author Marc A. Suchard
 */

public class CompleteHistorySimulatorParser extends AbstractXMLObjectParser {

    /* standard xml parser stuff follows */
    public static final String HISTORY_SIMULATOR = "completeHistorySimulator";
    public static final String TREE = "tree";
    public static final String REPLICATIONS = "replications";
    public static final String COUNTS = MarkovJumpsTreeLikelihoodParser.COUNTS;
    public static final String REWARDS = MarkovJumpsTreeLikelihoodParser.REWARDS;
    public static final String JUMP_TAG_NAME = MarkovJumpsTreeLikelihoodParser.JUMP_TAG_NAME;
    public static final String JUMP_TAG = MarkovJumpsTreeLikelihoodParser.JUMP_TAG;

    public static final String SYN_JUMPS = "reportSynonymousMutations";
    public static final String NON_SYN_JUMPS = "reportNonSynonymousMutations";
    public static final String SUM_SITES = "sumAcrossSites";

    public static final String BRANCH_SPECIFIC_SPECIFICATION = "branchSpecificSpecification";
    public static final String BRANCH_VARIABLE_PARAMETER = "variableParameter";
    public static final String VARIABLE_VALUE_PARAMETER = "valuesParameter";

    public static final String ANNOTATE_WITH_ALIGNMENT = "annotateWithAlignment";

    public static final String ALIGNMENT_ONLY = "alignmentOnly";
    
    public String getParserName() {
        return HISTORY_SIMULATOR;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        int nReplications = xo.getIntegerAttribute(REPLICATIONS);

        Tree tree = (Tree) xo.getChild(Tree.class);
        GammaSiteRateModel siteModel = (GammaSiteRateModel) xo.getChild(GammaSiteRateModel.class);
        BranchRateModel rateModel = (BranchRateModel) xo.getChild(BranchRateModel.class);

        if (rateModel == null)
            rateModel = new DefaultBranchRateModel();

        DataType dataType = siteModel.getSubstitutionModel().getDataType();

        String jumpTag = xo.getAttribute(JUMP_TAG_NAME, JUMP_TAG);

        boolean sumAcrossSites = xo.getAttribute(SUM_SITES, false);

        Parameter branchSpecificParameter = null;
        Parameter variableValueParameter = null;

        if (xo.hasChildNamed(BRANCH_SPECIFIC_SPECIFICATION)) {
            XMLObject cxo = xo.getChild(BRANCH_SPECIFIC_SPECIFICATION);
            branchSpecificParameter = (Parameter) cxo.getChild(BRANCH_VARIABLE_PARAMETER).getChild(Parameter.class);
            variableValueParameter = (Parameter) cxo.getChild(VARIABLE_VALUE_PARAMETER).getChild(Parameter.class);           
        }

        CompleteHistorySimulator history = new CompleteHistorySimulator(tree, siteModel, rateModel, nReplications,
                sumAcrossSites, branchSpecificParameter, variableValueParameter);

        XMLObject cxo = xo.getChild(COUNTS);
        if (cxo != null) {
            MarkovJumpsTreeLikelihoodParser.parseAllChildren(cxo, history, dataType.getStateCount(),
                    jumpTag, MarkovJumpsType.COUNTS, false);
        }

        cxo = xo.getChild(REWARDS);
        if (cxo != null) {
            MarkovJumpsTreeLikelihoodParser.parseAllChildren(cxo, history, dataType.getStateCount(),
                    jumpTag, MarkovJumpsType.REWARDS, false);
        }

        if (dataType instanceof Codons) {
            Codons codons = (Codons) dataType;
            if (xo.getAttribute(SYN_JUMPS, false)) {
                double[] synRegMatrix = CodonLabeling.getRegisterMatrix(CodonLabeling.SYN, codons, false); // use base 61
                Parameter registerParameter = new Parameter.Default(synRegMatrix);
                registerParameter.setId("S");
                history.addRegister(registerParameter, MarkovJumpsType.COUNTS, false);
            }
            if (xo.getAttribute(NON_SYN_JUMPS, false)) {
                double[] nonSynRegMatrix = CodonLabeling.getRegisterMatrix(CodonLabeling.NON_SYN, codons, false); // use base 61
                Parameter registerParameter = new Parameter.Default(nonSynRegMatrix);
                registerParameter.setId("N");
                history.addRegister(registerParameter, MarkovJumpsType.COUNTS, false);
            }
        }

        if (xo.getAttribute(ANNOTATE_WITH_ALIGNMENT,false)) {
            history.addAlignmentTrait();
        }

		boolean alignmentOnly = xo.getAttribute(ALIGNMENT_ONLY, false);
		if (dataType instanceof Codons && !alignmentOnly) {
			System.out.println("Codon models give exception when count statistics are done on them. "
							+ "You can supress this by setting alignmentOnly to true.");
		}

		if (alignmentOnly) {
			history.setAlignmentOnly();
		}
        
        history.simulate();
        return history;
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "A SequenceSimulator that generates random sequences for a given tree, siteratemodel and branch rate model";
    }

    public Class getReturnType() {
        return Alignment.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
            new ElementRule(Tree.class),
            new ElementRule(GammaSiteRateModel.class),
            new ElementRule(BranchRateModel.class, true),
            AttributeRule.newIntegerRule(REPLICATIONS),
            AttributeRule.newBooleanRule(SYN_JUMPS, true),
            AttributeRule.newBooleanRule(NON_SYN_JUMPS, true),
            AttributeRule.newBooleanRule(SUM_SITES, true),
            AttributeRule.newBooleanRule(ANNOTATE_WITH_ALIGNMENT, true),
            AttributeRule.newBooleanRule(ALIGNMENT_ONLY, true),
            new ElementRule(BRANCH_SPECIFIC_SPECIFICATION, new XMLSyntaxRule[] {
                    new ElementRule(VARIABLE_VALUE_PARAMETER, Parameter.class),
                    new ElementRule(BRANCH_VARIABLE_PARAMETER, Parameter.class),
            }, true),
    };
}