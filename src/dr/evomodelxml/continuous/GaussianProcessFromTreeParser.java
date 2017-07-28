/*
 * GaussianProcessFromTreeParser.java
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

package dr.evomodelxml.continuous;

import dr.evolution.tree.TreeTrait;
import dr.evomodel.continuous.FullyConjugateMultivariateTraitLikelihood;
import dr.evomodel.continuous.GaussianProcessFromTree;
import dr.evomodel.continuous.GibbsSampleMissingTraitsOperator;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.evomodel.treedatalikelihood.continuous.ContinuousDataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.continuous.TreeTipGaussianProcess;
import dr.math.distributions.GaussianProcessRandomGenerator;
import dr.xml.*;

/**
 * @author Marc A. Suchard
 * @author Max R. Tolkoff
 */
public class GaussianProcessFromTreeParser extends AbstractXMLObjectParser {

    public static final String GAUSSIAN_PROCESS_FROM_TREE = "gaussianProcessFromTree";
    public static final String MASK_TO_MISSING = "maskDrawToMissingOnly";

    private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
            new XORRule(
                    new ElementRule(FullyConjugateMultivariateTraitLikelihood.class),
                    new ElementRule(TreeDataLikelihood.class)
            ),
            AttributeRule.newBooleanRule(MASK_TO_MISSING, true),
    };

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        FullyConjugateMultivariateTraitLikelihood traitModel = (FullyConjugateMultivariateTraitLikelihood)
                xo.getChild(FullyConjugateMultivariateTraitLikelihood.class);

        if (traitModel != null) {
            return new GaussianProcessFromTree(traitModel);
        }

        TreeDataLikelihood treeDataLikelihood = (TreeDataLikelihood) xo.getChild(TreeDataLikelihood.class);
        ContinuousDataLikelihoodDelegate dataDelegate =
                GibbsSampleMissingTraitsOperator.parseContinuousDataLikelihoodDelegate(xo);

        boolean mask = xo.getAttribute(MASK_TO_MISSING, true);

        TreeTipGaussianProcess process = new TreeTipGaussianProcess(dataDelegate.getDataModel().getModelName(),
                treeDataLikelihood, dataDelegate, null, true);

        return process;

    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    @Override
    public String getParserDescription() {
        return "Returns a random draw of traits given a trait model and a prior";
    }

    @Override
    public Class getReturnType() {
        return GaussianProcessRandomGenerator.class;
    }

    @Override
    public String getParserName() {
        return GAUSSIAN_PROCESS_FROM_TREE;
    }
}
