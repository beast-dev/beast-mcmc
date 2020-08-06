/*
 * FullyConjugateTreeTipsPotentialDerivativeParser.java
 *
 * Copyright (c) 2002-2017 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.evomodelxml.continuous.hmc;

import dr.evomodel.continuous.MultivariateDiffusionModel;
import dr.evomodel.continuous.hmc.*;
import dr.inference.model.Parameter;
import dr.inferencexml.model.MaskedParameterParser;
import dr.xml.*;

/**
 * @author Aki Nishimura
 * @author Zhenyu Zhang
 * @author Marc A. Suchard
 */

public class TreePrecisionColumnProviderParser extends AbstractXMLObjectParser {

    private static final String PRODUCT_PROVIDER = "precisionColumnOnTree";
    private static final String KRONECKER_PRODUCT = "kroneckerProduct";
    private static final String EXTEND_TIP_BRANCH = "extendTipBranchTransformed";
    private static final String MASKING = MaskedParameterParser.MASKING;

    @Override
    public String getParserName() {
        return PRODUCT_PROVIDER;
    }

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        if (xo.hasChildNamed(KRONECKER_PRODUCT)) {

            XMLObject cxo = xo.getChild(KRONECKER_PRODUCT);

            TreePrecisionTraitProductProvider productProvider = (TreePrecisionTraitProductProvider)
                    cxo.getChild(TreePrecisionTraitProductProvider.class);

            if (productProvider.getDataModel().getTraitDimension() != 1) {
                throw new XMLParseException("Tree trait dimension should = 1 when used in a Kronecker product");
            }

            MultivariateDiffusionModel diffusionModel = (MultivariateDiffusionModel)
                    cxo.getChild(MultivariateDiffusionModel.class);
            boolean extendTipBranchTransformed = xo.getAttribute(EXTEND_TIP_BRANCH, false);

            return new TreeKroneckerPrecisionColumnProvider(productProvider, diffusionModel, extendTipBranchTransformed);

        } else {

            TreePrecisionTraitProductProvider productProvider = (TreePrecisionTraitProductProvider)
                    xo.getChild(TreePrecisionTraitProductProvider.class);

            return new TreePrecisionColumnProvider(productProvider);
        }
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new XORRule(
                    new ElementRule(TreePrecisionTraitProductProvider.class),
                    new ElementRule(KRONECKER_PRODUCT, new XMLSyntaxRule[] {
                            new ElementRule(TreePrecisionTraitProductProvider.class),
                            new ElementRule(MultivariateDiffusionModel.class),
                    })
            ),
            new ElementRule(MASKING,
                    new XMLSyntaxRule[]{
                            new ElementRule(Parameter.class)
                    }, true),
    };

    @Override
    public String getParserDescription() {
        return null;
    }

    @Override
    public Class getReturnType() {
        return TreePrecisionTraitProductProvider.class;
    }
}
