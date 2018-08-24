/*
 * OrderedLatentLiabilityTransformParser.java
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

package dr.evomodelxml.continuous;

import dr.evolution.datatype.DataType;
import dr.evolution.datatype.TwoStates;
import dr.evomodel.continuous.OrderedLatentLiabilityLikelihood;
import dr.inference.model.*;
import dr.inferencexml.model.MaskedParameterParser;
import dr.util.Transform;
import dr.xml.*;

import java.util.ArrayList;
import java.util.List;

/**
 * @author Marc A. Suchard
 */

public class OrderedLatentLiabilityTransformParser extends AbstractXMLObjectParser {

    public static final String NAME = "orderedLatentLiabilityTransform";

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        OrderedLatentLiabilityLikelihood likelihood = (OrderedLatentLiabilityLikelihood)
                xo.getChild(OrderedLatentLiabilityLikelihood.class);

        CompoundParameter parameter = likelihood.getTipTraitParameter();

        DataType dataType = likelihood.getPatternList().getDataType();

        if (!(dataType instanceof TwoStates)) {
            throw new XMLParseException("Liability transformation is currently only implemented for binary traits");
        }

        Parameter mask = null;
        if (xo.hasChildNamed(MaskedParameterParser.MASKING)) {
            mask = (Parameter) xo.getElementFirstChild(MaskedParameterParser.MASKING);
        }

        List<Transform> transforms = new ArrayList<Transform>();
        int index = 0;
        for (int tip = 0; tip < parameter.getParameterCount(); ++tip) {

            final int[] tipData = likelihood.getData(tip);

            for (int trait = 0; trait < tipData.length; ++trait) {
                int discreteState = tipData[trait];
                boolean valid = true;

                Transform transform;
                if (discreteState == 0) {
                    transform = Transform.LOG_NEGATE;
//                    transforms.add(Transform.LOG_NEGATE);

                    if (parameter.getParameterValue(index) >= 0.0) {
                        valid = false;
                    }
                } else if (discreteState == 1) {
                    transform = Transform.LOG;
//                    transforms.add(Transform.LOG);

                    if (parameter.getParameterValue(index) <= 0.0) {
                        valid = false;
                    }
                } else {
                    transform = Transform.NONE;
//                    transforms.add(Transform.NONE);
                }

                if (!valid) {
                    throw new XMLParseException("Incompatible binary trait and latent value in tip '" +
                            parameter.getParameter(tip).getId() + "'");
                }

                if (mask == null || mask.getParameterValue(index) == 1.0) {
                    transforms.add(transform);
                }

                ++index;
            }
        }

        return new Transform.Array(transforms, parameter);
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return new XMLSyntaxRule[] {
                new ElementRule(OrderedLatentLiabilityLikelihood.class),
                new ElementRule(MaskedParameterParser.MASKING,
                        new XMLSyntaxRule[]{
                                new ElementRule(Parameter.class)
                        }, true),
        };
    }

    @Override
    public String getParserDescription() {
        return null;
    }

    @Override
    public Class getReturnType() {
        return Transform.Collection.class;
    }

    @Override
    public String getParserName() {
        return NAME;
    }

}
