
/*
 * PrecisionColumnProviderParser.java
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
import dr.inference.hmc.PrecisionColumnProvider;
import dr.inference.model.MatrixParameterInterface;
import dr.xml.*;

/**
 * @author Zhenyu Zhang
 */

public class PrecisionColumnProviderParser extends AbstractXMLObjectParser {
    private static final String PRODUCT_PROVIDER = "precisionColumn";
    private static final String USE_CACHE = "useCache";

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        MatrixParameterInterface matrix = (MatrixParameterInterface) xo.getChild(MatrixParameterInterface.class);

        AutoRegressiveNormalDistributionModel ar = (AutoRegressiveNormalDistributionModel) xo.getChild(
                AutoRegressiveNormalDistributionModel.class);

        CompoundSymmetryNormalDistributionModel cs = (CompoundSymmetryNormalDistributionModel) xo.getChild(
                CompoundSymmetryNormalDistributionModel.class);

        boolean useCache = xo.getAttribute(USE_CACHE, true);

        if (matrix != null) {
            return new PrecisionColumnProvider.Generic(matrix, useCache);
        } else {
            if (ar != null) {
                return new PrecisionColumnProvider.AutoRegressive(ar, useCache);
            } else if (cs != null) {
                return new PrecisionColumnProvider.CompoundSymmetry(cs, useCache);
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
            AttributeRule.newBooleanRule(USE_CACHE, true),
    };

    @Override
    public String getParserDescription() {
        return null;
    }

    @Override
    public Class getReturnType() {
        return null;
    }

    @Override
    public String getParserName() {
        return PRODUCT_PROVIDER;
    }
}
