/*
 * CompoundSymmetryNormalDistributionModelParser.java
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

package dr.inferencexml.distribution;

import dr.inference.distribution.CompoundSymmetryNormalDistributionModel;
import dr.inference.model.Parameter;
import dr.xml.*;

/**
 * @author Zhenyu Zhang
 */

public class CompoundSymmetryNormalDistributionModelParser extends AbstractXMLObjectParser {

    public static final String NORMAL_DISTRIBUTION_MODEL = "compoundSymmetryNormalDistributionModel";
    private static final String DIMENSION = "dim";
    private static final String MARGINAL_VARIANCE = "variance";
    private static final String CORRELATION = "rho";

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        int dim = xo.getIntegerAttribute(DIMENSION);

        XMLObject cxo = xo.getChild(MARGINAL_VARIANCE);
        Parameter variance = (Parameter) cxo.getChild(Parameter.class);

        if (variance.getParameterValue(0) <= 0.0) {
            throw new XMLParseException("variance must be > 0.0");
        }

        cxo = xo.getChild(CORRELATION);
        Parameter rho = (Parameter) cxo.getChild(Parameter.class);

        if (Math.abs(rho.getParameterValue(0)) >= 1.0) {
            throw new XMLParseException("|Rho| must be < 1.0");
        }
        return new CompoundSymmetryNormalDistributionModel(dim, variance, rho);
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            AttributeRule.newIntegerRule(DIMENSION),
            new ElementRule(MARGINAL_VARIANCE,
                    new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
            new ElementRule(CORRELATION,
                    new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
    };

    @Override
    public String getParserDescription() {
        return null;
    }

    @Override
    public Class getReturnType() {
        return CompoundSymmetryNormalDistributionModel.class;
    }

    @Override
    public String getParserName() {
        return NORMAL_DISTRIBUTION_MODEL;
    }
}
