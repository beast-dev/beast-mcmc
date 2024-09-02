/*
 * PrecisionDataProductProviderParser.java
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

import dr.inference.distribution.AutoRegressiveNormalDistributionModel;
import dr.inference.distribution.CompoundSymmetryNormalDistributionModel;
import dr.inference.hmc.PrecisionMatrixVectorProductProvider;
import dr.inference.model.MatrixParameterInterface;
import dr.inference.model.Parameter;
import dr.inferencexml.model.MaskedParameterParser;
import dr.xml.*;

/**
 * @author Zhenyu Zhang
 */

public class PrecisionDataProductProviderParser extends AbstractXMLObjectParser {
    private static final String PRODUCT_PROVIDER = "precisionVectorProduct";
    private static final String MASKING = MaskedParameterParser.MASKING;
    private static final String TIME_GUESS = "roughTravelTimeGuess";
    private static final String THREAD_COUNT = "threadCount";

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        double roughTimeGuess = xo.getAttribute(TIME_GUESS, 1.0);

        MatrixParameterInterface matrix = (MatrixParameterInterface) xo.getChild(MatrixParameterInterface.class);

        AutoRegressiveNormalDistributionModel ar = (AutoRegressiveNormalDistributionModel) xo.getChild(
                AutoRegressiveNormalDistributionModel.class);

        CompoundSymmetryNormalDistributionModel cs = (CompoundSymmetryNormalDistributionModel) xo.getChild(
                CompoundSymmetryNormalDistributionModel.class);

        if (matrix != null) {
            return new PrecisionMatrixVectorProductProvider.Generic(matrix, roughTimeGuess);
        } else {
            if (ar != null) {
                return new PrecisionMatrixVectorProductProvider.AutoRegressive(ar, roughTimeGuess);
            } else if (cs != null) {
                return new PrecisionMatrixVectorProductProvider.CompoundSymmetry(cs, roughTimeGuess);
            } else {
                throw new RuntimeException("unrecognized type, must be ar or cs!");
            }
        }
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            AttributeRule.newIntegerRule(THREAD_COUNT, true),
            AttributeRule.newDoubleRule(TIME_GUESS, true),
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
        return PrecisionMatrixVectorProductProvider.class;
    }

    @Override
    public String getParserName() {
        return PRODUCT_PROVIDER;
    }
}
