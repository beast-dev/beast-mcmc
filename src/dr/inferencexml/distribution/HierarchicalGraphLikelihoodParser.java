/*
 * HierarchicalGraphLikelihoodParser.java
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

package dr.inferencexml.distribution;

import dr.inference.distribution.HierarchicalGraphLikelihood;
import dr.inference.model.Likelihood;
import dr.inference.model.MatrixParameter;
import dr.inference.model.Parameter;
import dr.xml.*;

/**
 * @author Gabriela Cybis
 */


public class HierarchicalGraphLikelihoodParser extends AbstractXMLObjectParser {

	

    public static final String HIERARCHICAL_INDICATOR = "hierarchicalIndicator";
    public static final String STRATA_INDICATOR = "strataIndicator";
    public static final String PROB = "prob";

    public String getParserName() {
        return HierarchicalGraphLikelihood.HIERARCHICAL_GRAPH_LIKELIHOOD;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        XMLObject cxo = xo.getChild(HIERARCHICAL_INDICATOR);
        Parameter hierarchicalIndicator = (Parameter) cxo.getChild(Parameter.class);

        cxo = xo.getChild(STRATA_INDICATOR);
        MatrixParameter strataIndicatorMatrix = new MatrixParameter("matrix");
        int dim = 0;
        for (int i = 0; i < cxo.getChildCount(); i++) {
            Parameter parameter = (Parameter) cxo.getChild(i);
            strataIndicatorMatrix.addParameter(parameter);
            if (i == 0)
                dim = parameter.getDimension();
            else if (dim != parameter.getDimension())
                throw new XMLParseException("All parameters must have the same dimension to construct a rectangular matrix");
        }
        
        if (hierarchicalIndicator.getDimension()!= strataIndicatorMatrix.getRowDimension())
        	throw new XMLParseException("Hierarchical and starta parameters don't have the same dimentions");
    
        
        cxo = xo.getChild(PROB);
        Parameter prob = (Parameter) cxo.getChild(Parameter.class);

        return new HierarchicalGraphLikelihood(hierarchicalIndicator, strataIndicatorMatrix, prob);

    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(HIERARCHICAL_INDICATOR,
                    new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
            new ElementRule(STRATA_INDICATOR,
                    new XMLSyntaxRule[]{new ElementRule(Parameter.class, 1, Integer.MAX_VALUE)}),
            new ElementRule(PROB,
                            new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
                    
    };

    public String getParserDescription() {
        return "Calculates the likelihood of strata graph given hierarchical graph and p.";
    }

    public Class getReturnType() {
        return Likelihood.class;
    }
}
