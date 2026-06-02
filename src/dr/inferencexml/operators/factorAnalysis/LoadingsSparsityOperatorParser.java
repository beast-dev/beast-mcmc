/*
 * LoadingsSparsityOperatorParser.java
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

package dr.inferencexml.operators.factorAnalysis;

import dr.inference.model.MatrixParameterInterface;
import dr.inference.operators.factorAnalysis.LoadingsGibbsTruncatedOperator;
import dr.inference.operators.factorAnalysis.LoadingsSparsityOperator;
import dr.xml.*;

public class LoadingsSparsityOperatorParser extends AbstractXMLObjectParser {
    public static final String LOADINGS_SPARSITY_OPERATOR = "loadingsSparsityOperator";
    public static final String WEIGHT = "weight";

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        LoadingsGibbsTruncatedOperator operator = (LoadingsGibbsTruncatedOperator) xo.getChild(LoadingsGibbsTruncatedOperator.class);
        MatrixParameterInterface sparseness = (MatrixParameterInterface) xo.getChild(MatrixParameterInterface.class);
        double weight = xo.getDoubleAttribute(WEIGHT);

        return new LoadingsSparsityOperator(weight, operator, sparseness);
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    XMLSyntaxRule[] rules = {
            new ElementRule(LoadingsGibbsTruncatedOperator.class),
            new ElementRule(MatrixParameterInterface.class),
    };

    @Override
    public String getParserDescription() {
        return "Sparseness operator for the sparseness on the loadings of a latent factor model";
    }

    @Override
    public Class getReturnType() {
        return LoadingsSparsityOperator.class;
    }

    @Override
    public String getParserName() {
        return LOADINGS_SPARSITY_OPERATOR;
    }
}
