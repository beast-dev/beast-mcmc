/*
 * LKJTransformParser.java
 *
 * Copyright (c) 2002-2018 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.util;

import dr.inference.model.Parameter;
import dr.xml.*;

import java.util.ArrayList;
import java.util.List;

/**
 * @author Paul Bastide
 */

public class LKJTransformParser extends AbstractXMLObjectParser {

    public static final String NAME = "LKJTransform";
    public static final String DIMENSION = "dimension";

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        int dim = xo.getIntegerAttribute(DIMENSION);
        int length = dim * (dim - 1) / 2;

        // Fisher Z  (unconstrained to constrained CPCs)
        List<Transform> transforms = new ArrayList<Transform>();
        Transform fisherZ = new Transform.Inverse(Transform.FISHER_Z);
        for (int i = 0; i < length; i++) {
            transforms.add(fisherZ);
        }
        Transform.Array fisherZTransforms = new Transform.Array(transforms, null);

        // LKJ inverse (constrained CPCs to correlation matrix)
        Transform.MultivariableTransform LKJTransform = new LKJTransformConstrained(dim);

        // Compose
        return new Transform.ComposeMultivariable(LKJTransform, fisherZTransforms);
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return new XMLSyntaxRule[]{
                AttributeRule.newIntegerRule(DIMENSION, true),
        };
    }

    @Override
    public String getParserDescription() {
        return "Returns a LKJ Transformation from unconstrained values to correlation matrix";
    }

    @Override
    public Class getReturnType() {
        return Transform.ComposeMultivariable.class;
    }

    @Override
    public String getParserName() {
        return NAME;
    }

}
