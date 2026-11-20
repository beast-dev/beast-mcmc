/*
 * ScaleStartingValueParser.java
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

package dr.inferencexml.model;

import dr.inference.model.Parameter;
import dr.xml.*;

import java.util.logging.Logger;

/**
 * Sets a scale parameter to the root mean square, sqrt(mean x^2), of every element of the other parameters, e.g. a
 * shrinkage prior's global scale from starting values of the parameters it shrinks. Place it after whatever sets
 * those starting values.
 * <p>
 * With an exponent a, the scale is instead (a * mean |x|^a)^(1/a), the mode of the scale tau given x under
 * independent exp(-|x / tau|^a) priors (a Bayesian bridge, ignoring its local scales and slab) and a flat prior on
 * tau^-a.
 * <pre>
 * &lt;scaleStartingValue&gt;
 *     &lt;scale&gt;&lt;parameter idref="globalScale"/&gt;&lt;/scale&gt;
 *     &lt;exponent&gt;&lt;parameter idref="exponent"/&gt;&lt;/exponent&gt;  (optional)
 *     &lt;parameter idref="rates"/&gt;
 * &lt;/scaleStartingValue&gt;
 * </pre>
 *
 * @author Marc A. Suchard
 */
public class ScaleStartingValueParser extends AbstractXMLObjectParser {

    private static final String PARSER_NAME = "scaleStartingValue";
    private static final String SCALE = "scale";
    private static final String EXPONENT = "exponent";

    public String getParserName() {
        return PARSER_NAME;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        Parameter scale = (Parameter) xo.getChild(SCALE).getChild(Parameter.class);

        double exponent = xo.hasChildNamed(EXPONENT) ?
                ((Parameter) xo.getChild(EXPONENT).getChild(Parameter.class)).getParameterValue(0) : Double.NaN;
        if (!Double.isNaN(exponent) && !(exponent > 0.0)) {
            throw new XMLParseException(EXPONENT + " must be positive");
        }
        double power = Double.isNaN(exponent) ? 2.0 : exponent;

        double sumOfPowers = 0.0;
        int count = 0;
        for (int i = 0; i < xo.getChildCount(); ++i) {
            if (xo.getChild(i) instanceof Parameter) {
                Parameter parameter = (Parameter) xo.getChild(i);
                for (int j = 0; j < parameter.getDimension(); ++j) {
                    sumOfPowers += Math.pow(Math.abs(parameter.getParameterValue(j)), power);
                }
                count += parameter.getDimension();
            }
        }

        double mean = sumOfPowers / count;
        double value = Double.isNaN(exponent) ? Math.sqrt(mean) : Math.pow(exponent * mean, 1.0 / exponent);
        if (!(value > 0.0) || Double.isInfinite(value)) {
            throw new XMLParseException("Cannot set " + scale.getId() + " to " + value);
        }

        for (int j = 0; j < scale.getDimension(); ++j) {
            scale.setParameterValueQuietly(j, value);
        }
        scale.fireParameterChangedEvent();

        Logger.getLogger("dr.inference").info("Set starting value of " + scale.getId() + " to " + value + ", " +
                (Double.isNaN(exponent) ? "the root mean square" : "the bridge scale mode (exponent " + exponent + ")") +
                " of " + count + " values");

        return value; // not the scale parameter: XMLParser would overwrite its id
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(SCALE, Parameter.class),
            new ElementRule(EXPONENT, Parameter.class, "bridge exponent; default: root mean square", true),
            new ElementRule(Parameter.class, 1, Integer.MAX_VALUE),
    };

    public String getParserDescription() {
        return "Sets a scale parameter to the root mean square, or the bridge scale mode, of other parameters";
    }

    public Class getReturnType() {
        return Double.class;
    }
}
