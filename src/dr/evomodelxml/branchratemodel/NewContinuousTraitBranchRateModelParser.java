/*
 * ContinuousTraitBranchRateModelParser.java
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

package dr.evomodelxml.branchratemodel;

import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.branchratemodel.NewContinuousTraitBranchRateModel;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.inference.model.Parameter;
import dr.xml.*;

import java.util.logging.Logger;

public class NewContinuousTraitBranchRateModelParser extends AbstractXMLObjectParser {

    public static final String TRAIT_BRANCH_RATES = "newContinuousTraitRateModel";
    public static final String TRAIT = "trait";
    public static final String DIMENSION = "dimension";
    public static final String RATE = BranchRateModel.RATE;
    public static final String RATIO = "ratio";

    public String getParserName() {
        return TRAIT_BRANCH_RATES;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        String trait = xo.getAttribute(TRAIT, "");
        int dimension = 0;
        if (xo.hasAttribute(DIMENSION)) {
            dimension = xo.getIntegerAttribute(DIMENSION) - 1;
        }

        TreeDataLikelihood traitLikelihood = (TreeDataLikelihood)
                xo.getChild(TreeDataLikelihood.class);
//        if (traitLikelihood != null)
//            trait = traitLikelihood.getTreeTraits()[dimension].getTraitName();

        Logger.getLogger("dr.evomodel").info("\nUsing trait '" + trait + "' as log rate estimates.");

        if (xo.hasChildNamed(RATE)) {
            Parameter rateParameter = (Parameter) xo.getElementFirstChild(RATE);
            Parameter ratioParameter = (Parameter) xo.getElementFirstChild(RATIO);

            return new NewContinuousTraitBranchRateModel(trait, rateParameter, ratioParameter);
        } else {

            return new NewContinuousTraitBranchRateModel(trait, traitLikelihood, dimension);
        }
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "This element returns an trait rate model." +
                "The branch rates are an average of the rates provided by a node trait.";
    }

    public Class getReturnType() {
        return NewContinuousTraitBranchRateModel.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
            AttributeRule.newStringRule(TRAIT, false, "The name of the trait that provides the log rates at nodes"),
            AttributeRule.newIntegerRule(DIMENSION, true, "The dimension that supplies the rate"),
            new ElementRule(RATE, Parameter.class, "The rate parameter", true),
            new ElementRule(RATIO, Parameter.class, "The ratio parameter", true),
            new ElementRule(TreeDataLikelihood.class, true),
    };


}
