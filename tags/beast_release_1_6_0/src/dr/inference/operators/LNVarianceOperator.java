/*
 * LNVarianceOperator.java
 *
 * Copyright (C) 2002-2006 Alexei Drummond and Andrew Rambaut
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

package dr.inference.operators;

import dr.inference.distribution.LogNormalDistributionModel;
import dr.inferencexml.operators.ScaleOperatorParser;
import dr.xml.*;

/**
 * This operator acts on the S parameter of the LogNormalDistribution but corrects the M parameter
 * so that the mean of the distribution is unchanged.
 *
 * @author Alexei Drummond
 * @version $Id: LNVarianceOperator.java,v 1.1 2004/12/16 14:17:21 alexei Exp $
 */
public class LNVarianceOperator extends ScaleOperator {

    public static final String LN_VARIANCE_OPERATOR = "LNVarianceOperator";


    private LogNormalDistributionModel lnd;

    public LNVarianceOperator(LogNormalDistributionModel lnd, double scaleFactor, double weight, CoercionMode mode) {

        //  super(lnd.getSParameter(), false, scaleFactor, weight, mode, null, 0.0);
        super(lnd.getSParameter(), false, 0, scaleFactor, mode, null, 0, false);
        this.lnd = lnd;
        setWeight(weight);
    }

    /**
     * Correct the M parameter of the lognormal distribution so that the mean is unchanged.
     */
    final void cleanupOperation(double newS, double oldS) {

        double newM = lnd.getM() + (oldS * oldS / 2.0) - (newS * newS / 2.0);
        lnd.setM(newM);
    }

    public static dr.xml.XMLObjectParser PARSER = new dr.xml.AbstractXMLObjectParser() {

        public String getParserName() {
            return LN_VARIANCE_OPERATOR;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            CoercionMode mode = CoercionMode.parseMode(xo);
            double weight = xo.getDoubleAttribute(WEIGHT);
            double scaleFactor = xo.getDoubleAttribute(ScaleOperatorParser.SCALE_FACTOR);

            if (scaleFactor <= 0.0 || scaleFactor >= 1.0) {
                throw new XMLParseException("scaleFactor must be between 0.0 and 1.0");
            }

            LogNormalDistributionModel lnd = (LogNormalDistributionModel) xo.getChild(LogNormalDistributionModel.class);

            return new LNVarianceOperator(lnd, scaleFactor, weight, mode);
        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public String getParserDescription() {
            return "This element returns a scale operator on a given parameter.";
        }

        public Class getReturnType() {
            return MCMCOperator.class;
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
                AttributeRule.newDoubleRule(ScaleOperatorParser.SCALE_FACTOR),
                AttributeRule.newDoubleRule(WEIGHT),
                AttributeRule.newBooleanRule(AUTO_OPTIMIZE, true),
                new ElementRule(LogNormalDistributionModel.class)
        };

    };


    public String toString() {
        return "LNVarianceOperator(" + lnd.getSParameter().getParameterName() + " [" + getScaleFactor() + ", " + (1.0 / getScaleFactor()) + "]";
    }
}
