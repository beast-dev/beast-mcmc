/*
 * GibbsSampleMissingTraitsOperator.java
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

package dr.evomodel.continuous;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeTrait;
import dr.evomodel.treedatalikelihood.DataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.evomodel.treedatalikelihood.continuous.ContinuousDataLikelihoodDelegate;
import dr.evomodelxml.treelikelihood.TreeTraitParserUtilities;
import dr.inference.model.Parameter;
import dr.inference.operators.MCMCOperator;
//import dr.inference.operators.OperatorFailedException;
import dr.inference.operators.SimpleMCMCOperator;
import dr.xml.*;

import java.util.ArrayList;
import java.util.List;

//import static dr.evomodel.treedatalikelihood.preorder.ConditionalOnPartiallyMissingTipsRealizedDelegate.PARTIAL;
import static dr.evomodelxml.treelikelihood.TreeTraitParserUtilities.MISSING;
import static dr.evomodelxml.treelikelihood.TreeTraitParserUtilities.TRAIT_NAME;

/**
 * @author Marc A. Suchard
 */
public class GibbsSampleMissingTraitsOperator extends SimpleMCMCOperator
//        implements GibbsOperator
{


    private static final String PARTIAL = "partial";

    final private TreeDataLikelihood treeLikelihood;
    final private TreeTrait treeTrait;
    final private Parameter parameter;
    final private Parameter missing;

    public GibbsSampleMissingTraitsOperator(TreeDataLikelihood treeLikelihood,
                                            TreeTrait treeTrait, Parameter parameter, Parameter missing) {
        super();

        this.treeLikelihood = treeLikelihood;
        this.treeTrait = treeTrait;

        this.parameter = parameter;
        this.missing = missing;
    }

    @Override
    public double doOperation() {//throws OperatorFailedException {

        final Tree tree = treeLikelihood.getTree();

        for (int i = 0; i < tree.getExternalNodeCount(); ++i) {
            final NodeRef node = tree.getExternalNode(i);

            treeTrait.getTrait(tree, node);
        }

        double[] tipValues = (double[]) treeTrait.getTrait(tree, null);
//        System.err.println(new dr.math.matrixAlgebra.Vector(tipValues));

        assert (tipValues.length == parameter.getDimension());

        for (int i = 0; i < parameter.getDimension(); ++i) {
            boolean fill = missing.getParameterValue(i) == 1;
            if (fill) {
                parameter.setParameterValue(i, tipValues[i]); // TODO Update silently?
            }
        }

        return 0; // TODO Not yet known ...
    }

    @Override
    public String getPerformanceSuggestion() {
        return "";
    }

    @Override
    public String getOperatorName() {
        return OPERATOR_NAME + "(" + "TODO" + ")";
    }

//    @Override
//    public int getStepCount() {
//        return 1;
//    }

    public static dr.xml.XMLObjectParser PARSER = new dr.xml.AbstractXMLObjectParser() {

        public String getParserName() {
            return OPERATOR_NAME;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            double weight = xo.getDoubleAttribute(MCMCOperator.WEIGHT);

            TreeDataLikelihood treeLikelihood = (TreeDataLikelihood) xo.getChild(TreeDataLikelihood.class);
//            ContinuousDataLikelihoodDelegate traitDelegate = parseContinuousDataLikelihoodDelegate(xo);

            TreeTrait treeTrait = parseTreeTrait(xo, PARTIAL);

            Parameter parameter = (Parameter) xo.getChild(Parameter.class);
            Parameter missing = (Parameter) xo.getElementFirstChild(MISSING);

            if (parameter.getDimension() != missing.getDimension()) {
                throw new XMLParseException("Unequal parameter lengths");
            }

            GibbsSampleMissingTraitsOperator operator = new GibbsSampleMissingTraitsOperator(treeLikelihood,
                    treeTrait,
                    parameter, missing);

            operator.setWeight(weight);

            return operator;

        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private final XMLSyntaxRule[] rules = {
                AttributeRule.newDoubleRule(MCMCOperator.WEIGHT),
                new ElementRule(TreeDataLikelihood.class),
                new ElementRule(Parameter.class),
                new ElementRule(TreeTraitParserUtilities.MISSING,
                        new XMLSyntaxRule[] {
                                new ElementRule(Parameter.class),
                        }),
                AttributeRule.newStringRule(TRAIT_NAME, true),
        };

        public String getParserDescription() {
            return "This element returns an independence coalescent sampler from a demographic model.";
        }

        public Class getReturnType() {
            return GibbsSampleMissingTraitsOperator.class;
        }

    };

    public static ContinuousDataLikelihoodDelegate parseContinuousDataLikelihoodDelegate(XMLObject xo)
            throws XMLParseException {

        TreeDataLikelihood treeLikelihood = (TreeDataLikelihood) xo.getChild(TreeDataLikelihood.class);
        DataLikelihoodDelegate delegate = treeLikelihood.getDataLikelihoodDelegate();

        if (!(delegate instanceof ContinuousDataLikelihoodDelegate)) {
            throw new XMLParseException("Only implemented for multivariate trait diffusion models");
        }

        return (ContinuousDataLikelihoodDelegate) delegate;
    }

    private static TreeTrait parseTreeTrait(XMLObject xo, String prefix) throws XMLParseException {

        TreeDataLikelihood treeLikelihood = (TreeDataLikelihood) xo.getChild(TreeDataLikelihood.class);

        String traitName = null;
        if (xo.hasAttribute(TRAIT_NAME)) {
            traitName = (String) xo.getAttribute(TRAIT_NAME);
        } else {
            TreeTrait[] traits = treeLikelihood.getTreeTraits();
            List<String> traitNames = matchedTraitNames(traits, prefix + ".");

            if (traitNames.size() == 1) {
                traitName = traitNames.get(0);
            }
        }

        TreeTrait treeTrait = treeLikelihood.getTreeTrait(traitName);

        if (treeTrait == null) {
            throw new XMLParseException("Unknown partially observed tree trait");
        }

        return treeTrait;
    }

    private static List<String> matchedTraitNames(final TreeTrait[] traits, final String prefix) {
        List<String> names = new ArrayList<String>();
        for (TreeTrait trait : traits) {
            if (trait.getTraitName().startsWith(prefix)) {
                names.add(trait.getTraitName());
            }
        }
        return names;
    }

    public static final String OPERATOR_NAME = "gibbsSampleMissingTraitsOperator";

}
