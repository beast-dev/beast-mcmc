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

import dr.xml.*;

import java.util.ArrayList;
import java.util.List;

/**
 * @author Paul Bastide
 */

public class LKJTransformParser extends AbstractXMLObjectParser {

    public static final String NAME = "LKJTransform";
    public static final String DIMENSION = "dimension";
    private static final String CHOLESKY = "cholesky";
    public static final String WITH_DIAGONALS = "withDiagonals";

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        int dim = xo.getIntegerAttribute(DIMENSION);
        int length = dim * (dim - 1) / 2;

        // Fisher Z  (constrained CPCs to unconstrained)
        List<Transform> transforms = new ArrayList<Transform>();
//        Transform fisherZ = new Transform.Inverse(Transform.FISHER_Z);
        for (int i = 0; i < length; i++) {
            transforms.add(Transform.FISHER_Z);
        }
        Transform.Array fisherZTransforms = new Transform.Array(transforms, null);

        // LKJ (constrained CPCs to (cholesky of) correlation matrix)
        Transform.MultivariableTransform LKJTransform;
        boolean cholesky = xo.getAttribute(CHOLESKY, true);
        if (cholesky) {
            LKJTransform = new LKJCholeskyTransformConstrained(dim);
        } else {
            LKJTransform = new LKJTransformConstrained(dim);
// Should work too, but a bit slower:
//            LKJTransform = new Transform.ComposeMultivariable(
//                    new LKJCholeskyTransformConstrained(dim),
//                    new CorrelationToCholesky(dim));
        }

        // Compose
        // return new Transform.ComposeMultivariable(fisherZTransforms, LKJTransform);
        boolean withDiag = xo.getAttribute(WITH_DIAGONALS, false);

        if (!withDiag) {
            return new Transform.ComposeMultivariable(fisherZTransforms, LKJTransform);
        } else {

            LKJCholeskyTransformConstrainedWithDiag LKJwithLOG;

            if (cholesky) {
                LKJwithLOG = new LKJCholeskyTransformConstrainedWithDiag(dim);
            } else {
                throw new RuntimeException("Not yet implemented");
            }

            for (int i = 0; i < dim; i++) {
                transforms.add(Transform.LOG);
            }

            Transform.Array fisherZTransformsWithNull = new Transform.Array(transforms, null);
            Transform jointTrans = new Transform.ComposeMultivariable(fisherZTransformsWithNull, LKJwithLOG);

            return jointTrans;
        }
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return new XMLSyntaxRule[]{
                AttributeRule.newIntegerRule(DIMENSION, false),
                AttributeRule.newBooleanRule(CHOLESKY, true),
        };
    }

    @Override
    public String getParserDescription() {
        return "Returns a LKJ Transformation from correlation matrix to unconstrained values";
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
