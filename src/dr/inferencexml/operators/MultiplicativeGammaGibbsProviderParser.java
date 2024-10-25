/*
 * MultiplicativeGammaGibbsProviderParser.java
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

package dr.inferencexml.operators;

import dr.evomodel.continuous.MatrixShrinkageLikelihood;
import dr.inference.model.CompoundParameter;
import dr.inference.model.MatrixParameterInterface;
import dr.inference.model.Parameter;
import dr.inference.operators.repeatedMeasures.GammaGibbsProvider;
import dr.inference.operators.repeatedMeasures.MultiplicativeGammaGibbsHelper;
import dr.xml.*;

import java.util.ArrayList;
import java.util.List;

/**
 * @author Gabriel Hassler
 * @author Marc A. Suchard
 */

public class MultiplicativeGammaGibbsProviderParser extends AbstractXMLObjectParser {

    private static final String MULTIPLICATIVE_PROVIDER = "multiplicativeGammaGibbsProvider";


    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        Parameter multParam = (Parameter) xo.getChild(Parameter.class);

        MultiplicativeGammaGibbsHelper helper =
                (MultiplicativeGammaGibbsHelper) xo.getChild(MultiplicativeGammaGibbsHelper.class);


        int k = multParam.getDimension();

        if (helper.getColumnDimension() != k) {
            throw new XMLParseException("Dimension mismatch: the parameter with id `" + multParam.getId() + "` has dimension " + k +
                    ", while the helper has " + helper.getColumnDimension() + " columns.");
        }


        return new GammaGibbsProvider.MultiplicativeGammaGibbsProvider(multParam, helper);
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return new XMLSyntaxRule[]{
                new ElementRule(Parameter.class),
                new ElementRule(MultiplicativeGammaGibbsHelper.class),
        };
    }

    @Override
    public String getParserDescription() {
        return "Returns sufficient statistics for multiplicative gamma distribution.";
    }

    @Override
    public Class getReturnType() {
        return GammaGibbsProvider.MultiplicativeGammaGibbsProvider.class;
    }

    @Override
    public String getParserName() {
        return MULTIPLICATIVE_PROVIDER;
    }
}
