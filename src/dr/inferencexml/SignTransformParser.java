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

package dr.inferencexml;

import dr.inference.model.Bounds;
import dr.inference.model.Parameter;
import dr.util.Transform;
import dr.xml.*;

import java.util.ArrayList;
import java.util.List;

/**
 * @author Marc A. Suchard
 */

public class SignTransformParser extends AbstractXMLObjectParser {

    public static final String NAME = "signTransform";

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        Parameter parameter = (Parameter) xo.getChild(Parameter.class);
        Bounds<Double> bounds = parameter.getBounds();

        List<Transform> transforms = new ArrayList<Transform>();

        for (int i = 0; i < parameter.getDimension(); i++) { // TODO much better checking is necessary (here we assumed bounds <0 or >0 )
            if (bounds.getLowerLimit(i) == 0.0) {
                transforms.add(Transform.LOG);
            } else if (bounds.getUpperLimit(i) == 0.0) {
                transforms.add(Transform.LOG_NEGATE);
            } else {
                transforms.add(Transform.NONE);
            }
//            transforms.add(Transform.LOG); // TODO What is this doing here?  REMOVE?
        }

//        for (int i = 0; i < parameter.getDimension(); i++) {
//            Transform transform;
//            transform = Transform.LOG;
//            transforms.add(transform);
//        }

        return new Transform.Array(transforms, parameter);
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return new XMLSyntaxRule[] {
                new ElementRule(Parameter.class),

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
