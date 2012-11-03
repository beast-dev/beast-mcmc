/*
 * MultivariateNormalDistributionModelParser.java
 *
 * Copyright (c) 2002-2012 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

import dr.inference.distribution.MultivariateDistributionLikelihood;
import dr.inference.distribution.MultivariateNormalDistributionModel;
import dr.inference.model.MatrixParameter;
import dr.inference.model.Parameter;
import dr.xml.*;

/**
 * Reads a normal distribution model from a DOM Document element.
 */
public class MultivariateNormalDistributionModelParser extends AbstractXMLObjectParser {

    public static final String NORMAL_DISTRIBUTION_MODEL = "multivariateNormalDistributionModel";
//    public static final String MEAN = "mean";
//    public static final String PREC = "precision";

    public String getParserName() {
        return NORMAL_DISTRIBUTION_MODEL;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        XMLObject cxo = xo.getChild(MultivariateDistributionLikelihood.MVN_MEAN);
        Parameter mean = (Parameter) cxo.getChild(Parameter.class);

        cxo = xo.getChild(MultivariateDistributionLikelihood.MVN_PRECISION);
        MatrixParameter precision = (MatrixParameter) cxo.getChild(MatrixParameter.class);

        if (mean.getDimension() != precision.getRowDimension() ||
                mean.getDimension() != precision.getColumnDimension())
            throw new XMLParseException("Mean and precision have wrong dimensions in " + xo.getName() + " element");

        return new MultivariateNormalDistributionModel(mean, precision);
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(MultivariateDistributionLikelihood.MVN_MEAN,
                    new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
            new ElementRule(MultivariateDistributionLikelihood.MVN_PRECISION,
                    new XMLSyntaxRule[]{new ElementRule(MatrixParameter.class)}),
//            new ElementRule(MEAN,
//                    new XMLSyntaxRule[]{
//                            new XORRule(
//                                    new ElementRule(Parameter.class),
//                                    new ElementRule(Double.class)
//                            )}
//            ),
//            new XORRule(
//                    new ElementRule(STDEV,
//                            new XMLSyntaxRule[]{
//                                    new XORRule(
//                                            new ElementRule(Parameter.class),
//                                            new ElementRule(Double.class)
//                                    )}
//                    ),
//                    new ElementRule(PREC,
//                            new XMLSyntaxRule[]{
//                                    new XORRule(
//                                            new ElementRule(Parameter.class),
//                                            new ElementRule(Double.class)
//                                    )}
//                    )
//            )
    };

    public String getParserDescription() {
        return "Describes a normal distribution with a given mean and standard deviation " +
                "that can be used in a distributionLikelihood element";
    }

    public Class getReturnType() {
        return MultivariateNormalDistributionModel.class;
    }

}
