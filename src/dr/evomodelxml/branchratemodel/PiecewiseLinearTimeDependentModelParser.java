/*
 * PiecewiseLinearTimeDependentModelParser.java
 *
 * Copyright (c) 2002-2022 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.evomodelxml.branchratemodel;

import dr.evomodel.branchratemodel.PiecewiseLinearTimeDependentModel;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.Parameter;
import dr.xml.*;

import java.util.logging.Logger;

/**
 * @author Marc Suchard
 * @author Philippe Lemey
 */
public class PiecewiseLinearTimeDependentModelParser extends AbstractXMLObjectParser {

    public static final String EFFECT_NAME = "piecewiseLinearTimeEffect";
    private static final String RATES = "rates";
    private static final String TRANSITION_TIME = "epochStart";
    private static final String EPOCH_LENGTH = "epochLength";
    private static final String SCALE = "scale";

    public String getParserName() {
        return EFFECT_NAME;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        Logger.getLogger("dr.evomodel").info("\nUsing a piece-wise linear rate model.");

        XMLObject cxo = xo.getChild(RATES);
        Parameter historicalRate = (Parameter) cxo.getChild(0);
        Parameter currentRate = (Parameter) cxo.getChild(1);

        Parameter firstTransitionTime = (Parameter) xo.getElementFirstChild(TRANSITION_TIME);
        Parameter epochLength = (Parameter) xo.getElementFirstChild(EPOCH_LENGTH);

        TreeModel tree = (TreeModel) xo.getChild(TreeModel.class);

        PiecewiseLinearTimeDependentModel.Scale scale;

        try {
            scale = PiecewiseLinearTimeDependentModel.Scale.parse((String) xo.getAttribute(SCALE));
            if (scale == null) {
                throw new Exception();
            }
        } catch (Exception e) {
            throw new XMLParseException("Unknown scale");
        }

        return new PiecewiseLinearTimeDependentModel(tree,
                new PiecewiseLinearTimeDependentModel.ParameterPack(
                        historicalRate, currentRate,
                        firstTransitionTime, epochLength),
                scale);
    }
    
    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "A parser";
    }

    public Class getReturnType() {
        return PiecewiseLinearTimeDependentModel.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(TreeModel.class),
            new ElementRule(RATES, new XMLSyntaxRule[] {
                    new ElementRule(Parameter.class, 2, 2),
            }),
            new ElementRule(TRANSITION_TIME, Parameter.class),
            new ElementRule(EPOCH_LENGTH, Parameter.class),
            AttributeRule.newStringRule(SCALE),
    };
}
