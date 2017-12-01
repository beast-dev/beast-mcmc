/*
 * MaskedGradientParser.java
 *
 * Copyright (c) 2002-2017 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.inferencexml.hmc;

import dr.inference.model.Bounds;
import dr.inference.model.Parameter;
import dr.math.matrixAlgebra.Vector;
import dr.xml.*;
import org.ejml.alg.dense.decomposition.eig.watched.WatchedDoubleStepQREigenvector;

import java.util.Arrays;

/**
 * @author Marc A. Suchard
 */
public class RotationTranslationMaskParser extends AbstractXMLObjectParser {

    private final static String MASK = "rotationalTranslationalMask";
    private final static String DIMENSION = "dimension";
    private final static String RESET = "reset";

    @Override
    public String getParserName() {
        return MASK;
    }

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        Parameter parameter = (Parameter) xo.getChild(Parameter.class);
        int dim = xo.getIntegerAttribute(DIMENSION);

        if (parameter.getDimension() % dim != 0) {
            throw new XMLParseException("Dimension and parameter length are not divisible");
        }

        boolean reset = xo.getAttribute(RESET, true);

        Parameter mask = new Parameter.Default(parameter.getDimension(), 1.0);

        int offset = 0;

        // Translational invariance
        for (int i = 0; i < dim; ++i) {
            mask.setParameterValue(offset, 0.0);
            if (reset) {
                parameter.setParameterValue(offset, 0.0);
            }
            ++offset;
        }

        // Reflection invariance
        if (reset) {
            parameter.setParameterValue(offset,
                    Math.abs(parameter.getParameterValue(offset)));
        }

        // TODO The following does not yet work for masked parameter
//        double[] lower = new double[parameter.getDimension()];
//        Arrays.fill(lower, Double.NEGATIVE_INFINITY);
//
//        double[] upper = new double[parameter.getDimension()];
//        Arrays.fill(upper, Double.POSITIVE_INFINITY);
//
//        lower[offset] = 0.0;
//        parameter.addBounds(new Parameter.DefaultBounds(upper, lower));

        // Rotational invariance
        for (int column = 1; column < dim; ++column) {

            for (int i = 0; i < dim; ++i) {
                if (i >= column) {
                    mask.setParameterValue(offset, 0.0);
                    if (reset) {
                        parameter.setParameterValue(offset, 0.0);
                    }
                }
                ++offset;
            }
        }
        
        return mask;
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            AttributeRule.newIntegerRule(DIMENSION),
            AttributeRule.newBooleanRule(RESET, true),
            new ElementRule(Parameter.class),
    };

    @Override
    public String getParserDescription() {
        return null;
    }

    @Override
    public Class getReturnType() {
        return Parameter.class;
    }
}
