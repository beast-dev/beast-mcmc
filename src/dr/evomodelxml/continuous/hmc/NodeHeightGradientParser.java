/*
 * NodeHeightGradientParser.java
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

package dr.evomodelxml.continuous.hmc;

import dr.evolution.tree.Tree;
import dr.evomodel.branchratemodel.ArbitraryBranchRates;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.treedatalikelihood.BeagleDataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.DataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.evomodel.treedatalikelihood.continuous.ContinuousDataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.discrete.NodeHeightGradientForDiscreteTrait;
import dr.evomodelxml.treelikelihood.TreeTraitParserUtilities;
import dr.inference.hmc.GradientWrtParameterProvider;
import dr.inference.hmc.JointGradient;
import dr.inference.model.CompoundLikelihood;
import dr.inference.model.Likelihood;
import dr.inference.model.Parameter;
import dr.xml.*;

import java.util.ArrayList;
import java.util.List;

import static dr.evomodelxml.treelikelihood.TreeTraitParserUtilities.DEFAULT_TRAIT_NAME;

/**
 * @author Marc A. Suchard
 * @author Xiang Ji
 */

public class NodeHeightGradientParser extends AbstractXMLObjectParser {

    private static final String NAME = "nodeHeightGradient";
    private static final String TRAIT_NAME = TreeTraitParserUtilities.TRAIT_NAME;

    public String getParserName(){ return NAME; }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        String traitName = xo.getAttribute(TRAIT_NAME, DEFAULT_TRAIT_NAME);

        final TreeDataLikelihood treeDataLikelihood = (TreeDataLikelihood) xo.getChild(TreeDataLikelihood.class);
        if (treeDataLikelihood != null) {
            return parseTreeDataLikelihood(treeDataLikelihood, traitName);
        }

        // several partitions on one tree: sum their node-height gradients
        CompoundLikelihood compoundLikelihood = (CompoundLikelihood) xo.getChild(CompoundLikelihood.class);
        List<GradientWrtParameterProvider> gradients = new ArrayList<>();
        Tree tree = null;
        for (Likelihood likelihood : compoundLikelihood.getLikelihoods()) {
            if (!(likelihood instanceof TreeDataLikelihood)) {
                throw new XMLParseException("Node-height gradients need tree data likelihoods; " +
                        likelihood.getId() + " is not one");
            }
            TreeDataLikelihood partition = (TreeDataLikelihood) likelihood;
            if (tree == null) {
                tree = partition.getTree();
            } else if (partition.getTree() != tree) {
                throw new XMLParseException("Node-height gradients over a compound likelihood need a single shared " +
                        "tree; " + partition.getId() + " uses another");
            }
            gradients.add(parseTreeDataLikelihood(partition, traitName));
        }
        if (gradients.isEmpty()) {
            throw new XMLParseException("The compound likelihood " + compoundLikelihood.getId() + " is empty");
        }
        return gradients.size() == 1 ? gradients.get(0) : new JointGradient(gradients);
    }

    private static GradientWrtParameterProvider parseTreeDataLikelihood(TreeDataLikelihood treeDataLikelihood,
                                                                        String traitName) throws XMLParseException {
        BranchRateModel branchRateModel = treeDataLikelihood.getBranchRateModel();

        Parameter branchRates = null;
        if (branchRateModel instanceof ArbitraryBranchRates) {
            branchRates = ((ArbitraryBranchRates) branchRateModel).getRateParameter();
        }

        DataLikelihoodDelegate delegate = treeDataLikelihood.getDataLikelihoodDelegate();

        if (delegate instanceof ContinuousDataLikelihoodDelegate) {

//                ContinuousDataLikelihoodDelegate continuousData = (ContinuousDataLikelihoodDelegate) delegate;
//                return new NodeHeightGradient(traitName, treeDataLikelihood, continuousData, branchRates);
            throw new XMLParseException("Not yet implemented! ");

        } else if (delegate instanceof BeagleDataLikelihoodDelegate) {

            BeagleDataLikelihoodDelegate beagleData = (BeagleDataLikelihoodDelegate) delegate;
            if (!beagleData.isUsePreOrder()) {
                throw new XMLParseException("Node-height gradients need pre-order traversals: give " +
                        treeDataLikelihood.getId() + " the attribute usePreOrder=\"true\"");
            }
            return new NodeHeightGradientForDiscreteTrait(traitName, treeDataLikelihood, beagleData, branchRates);

        } else {
            throw new XMLParseException("Unknown likelihood delegate type");
        }
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            AttributeRule.newStringRule(TRAIT_NAME),
            new XORRule(
                    new ElementRule(TreeDataLikelihood.class),
                    new ElementRule(CompoundLikelihood.class)
            ),
    };

    @Override
    public String getParserDescription() {
        return null;
    }

    @Override
    public Class getReturnType() {
        return GradientWrtParameterProvider.class;
    }
}
