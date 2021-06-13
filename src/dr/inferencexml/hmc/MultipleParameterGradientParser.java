/*
 * BranchSpecificOptimaGradientParser.java
 *
 * Copyright (c) 2002-2019 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.inferencexml.hmc;

import dr.evomodel.branchratemodel.ArbitraryBranchRates;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.evomodel.treedatalikelihood.hmc.AbstractDiffusionGradient;
import dr.inference.hmc.GradientWrtParameterProvider;
import dr.inference.hmc.JointGradient;
import dr.inference.hmc.MultipleParameterGradient;
import dr.inference.model.CompoundParameter;
import dr.xml.*;

import java.util.ArrayList;
import java.util.List;

/**
 * @author Alexander Fisher
 */

public class MultipleParameterGradientParser extends AbstractXMLObjectParser {

    private final static String MULTI_PARAMETER_GRADIENT = "multipleParameterGradient";

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        CompoundParameter compoundParameter = new CompoundParameter(null);

        List<GradientWrtParameterProvider> gradients = xo.getAllChildren(GradientWrtParameterProvider.class);
        List<GradientWrtParameterProvider> derivativeList = new ArrayList<>();

        if (gradients != null) {
            for (GradientWrtParameterProvider grad : gradients) {
//                derivationParametersList.add(grad.getDerivationParameter());
                compoundParameter.addParameter(grad.getParameter());
                derivativeList.add(grad);
            }
        }

//        for (int i = 0; i < xo.getChildCount(); i++) {
//            GradientWrtParameterProvider grad = (GradientWrtParameterProvider) xo.getChild(i);
//            derivativeList.add(grad);
//        }

        return new MultipleParameterGradient(derivativeList, compoundParameter);
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
//            new ElementRule(TreeDataLikelihood.class),
//            new ElementRule(ArbitraryBranchRates.class, 1, Integer.MAX_VALUE),
            new ElementRule(GradientWrtParameterProvider.class, 1, Integer.MAX_VALUE),
    };
//
//    @Override
//    public XMLSyntaxRule[] getSyntaxRules() {
//        return rules;
//    }
//
//    private final XMLSyntaxRule[] rules = {
//            new ElementRule(AbstractDiffusionGradient.class, 1, Integer.MAX_VALUE),
//    };

    @Override
    public String getParserDescription() {
        return null;
    }

    @Override
    public Class getReturnType() {
        return MultipleParameterGradient.class;
    }

    @Override
    public String getParserName() {
        return MULTI_PARAMETER_GRADIENT;
    }
}
