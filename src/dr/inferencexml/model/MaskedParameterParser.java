/*
 * MaskedParameterParser.java
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

package dr.inferencexml.model;

import dr.inference.model.MaskedParameter;
import dr.inference.model.Parameter;
import dr.xml.*;

import java.util.logging.Logger;

/**
 * @author Marc A. Suchard
 */
public class MaskedParameterParser extends AbstractXMLObjectParser {

    public static final String MASKED_PARAMETER = "maskedParameter";
    public static final String MASKING = "mask";
    public static final String COMPLEMENT = "complement";
    public static final String FROM = "from";
    public static final String TO = "to";
    public static final String EVERY = "every";
    public static final String BUILD = "build";
    private static final String SIGNAL_DEPENDENTS = "signalDependents";
    private static final String IS_NA_MISSING = "isNaMissing";

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        Parameter parameter = (Parameter) xo.getChild(Parameter.class);
        Parameter mask;

        XMLObject cxo = xo.getChild(MASKING);
        if (cxo != null) {
            mask = (Parameter) cxo.getChild(Parameter.class);
        } else if (xo.getAttribute(BUILD, false)) {
            boolean isNaMissing = xo.getAttribute(IS_NA_MISSING, false);
            mask = new Parameter.Default(parameter.getDimension(), 0.0);
            for (int i = 0; i < parameter.getDimension(); i++) {
                if (isNaMissing && Double.isNaN(parameter.getParameterValue(i))) {
                    mask.setParameterValue(i, 1.0);
                    parameter.setParameterValue(i, 0.0);
                    Logger.getLogger("dr.inferencexml.model").info("Setting dim " + (i + 1) + " in " +
                            parameter.getId() + " to 0.0");

                }
                if (!isNaMissing && parameter.getParameterValue(i) == 0) {
                    mask.setParameterValue(i, 1.0);
                }
            }
        }
        else {
            int from = xo.getAttribute(FROM, 1) - 1;
            int to = xo.getAttribute(TO, parameter.getDimension()) - 1;
            int every = xo.getAttribute(EVERY, 1);

            if (from < 0) throw new XMLParseException("illegal 'from' attribute in " + MASKED_PARAMETER
                    + " element");
            if (to < 0 || to < from) throw new XMLParseException("illegal 'to' attribute in "
                    + MASKED_PARAMETER + " element");
            if (every <= 0) throw new XMLParseException("illegal 'every' attribute in " + MASKED_PARAMETER
                    + " element");

            mask = new Parameter.Default(parameter.getDimension(), 0.0);
            for (int i = from; i <= to; i += every) {
                mask.setParameterValue(i, 1.0);
            }
        }

        if (mask.getDimension() != parameter.getDimension())
            throw new XMLParseException("dim(" + parameter.getId() + ":" + parameter.getDimension()
                    + ") != dim(" + mask.getId() + ":" + mask.getDimension() + ")");

        boolean signal = xo.getAttribute(SIGNAL_DEPENDENTS, true);

        MaskedParameter maskedParameter = new MaskedParameter(parameter,
                signal ? MaskedParameter.Signaling.NORMAL : MaskedParameter.Signaling.NO_DEPENDENT);
        boolean ones = !xo.getAttribute(COMPLEMENT, false);
        maskedParameter.addMask(mask, ones);

        return maskedParameter;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(Parameter.class),
            new ElementRule(MASKING,
                    new XMLSyntaxRule[]{
                            new ElementRule(Parameter.class)
                    }, true),
            AttributeRule.newBooleanRule(COMPLEMENT, true),
            AttributeRule.newBooleanRule(BUILD, true),
            AttributeRule.newIntegerRule(FROM, true),
            AttributeRule.newIntegerRule(TO, true),
            AttributeRule.newIntegerRule(EVERY, true),
            AttributeRule.newBooleanRule(SIGNAL_DEPENDENTS, true),
            AttributeRule.newBooleanRule(IS_NA_MISSING, true),
    };

    public String getParserDescription() {
        return "A masked parameter.";
    }

    public Class getReturnType() {
        return Parameter.class;
    }

    public String getParserName() {
        return MASKED_PARAMETER;
    }
}
