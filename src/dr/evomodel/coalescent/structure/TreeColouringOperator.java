/*
 * TreeColouringOperator.java
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

package dr.evomodel.coalescent.structure;

import dr.evolution.colouring.ColourSampler;
import dr.inference.operators.SimpleMCMCOperator;
import dr.xml.*;

/**
 * @author Alexei Drummond
 *
 * @version $Id: TreeColouringOperator.java,v 1.8 2006/08/12 12:55:44 gerton Exp $
 */
public class TreeColouringOperator extends SimpleMCMCOperator {

    public static final String TREE_COLOURING_OPERATOR = "treeColouringOperator";

    ColourSamplerModel colouringModel;
    ColourSampler colourSampler;

    public TreeColouringOperator(ColourSamplerModel colouringModel) {

        this.colouringModel = colouringModel;
    }

    /**
     * @return a short descriptive message of the performance of this operator.
     */
    public String getPerformanceSuggestion() {
        return "This operator cannot be optimized";
    }

    public String getOperatorName() {
        return "twoColourTree";
    }

    /**
     * Called by operate(), does the actual operation.
     *
     * @return the hastings ratio
     */
    public double doOperation() {

        double logP = colouringModel.getTreeColouringWithProbability().getLogProbabilityDensity();

        colouringModel.resample();

        double logQ = colouringModel.getTreeColouringWithProbability().getLogProbabilityDensity();

        return logP - logQ;
    }

    /**
     * Reads a distribution likelihood from a DOM Document element.
     */
    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public String getParserName() { return TREE_COLOURING_OPERATOR; }

        public Object parseXMLObject(XMLObject xo) {

            ColourSamplerModel colourSamplerModel = (ColourSamplerModel)xo.getChild(ColourSamplerModel.class);

            return new TreeColouringOperator(colourSamplerModel);
        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public XMLSyntaxRule[] getSyntaxRules() { return rules; }

        private XMLSyntaxRule[] rules = new XMLSyntaxRule[] {
            new ElementRule(ColourSamplerModel.class),
        };

        public String getParserDescription() {
            return "A tree colouring model.";
        }

        public Class getReturnType() { return TreeColouringOperator.class; }
    };


}
