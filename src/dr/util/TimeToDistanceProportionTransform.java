/*
 * TimeToDistanceProportionTransform.java
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

package dr.util;

import dr.xml.AbstractXMLObjectParser;
import dr.xml.XMLObject;
import dr.xml.XMLParseException;
import dr.xml.XMLSyntaxRule;

public class TimeToDistanceProportionTransform extends Transform.MultivariateTransform {
    public static String NAME = "TimeToDistanceProportionTransform";

    // The variables are assumed to be in the following order:
    // 0: the proportion (in time for forward transform, in distance for reverse transform)
    // 1: the log-scale rate in the ancestral-model portion of the branch
    // 2: the log-scale rate in the descendant-model portion of the branch
    // The input is these as a length-3 vector, the output is a length-1 vector of the transformed proportion
    TimeToDistanceProportionTransform() {
        super(3, 1);
    }

    @Override
    public double[] inverse(double[] values, int from, int to, double sum) {
        throw new RuntimeException("Not yet implemented");
    }

    @Override
    public double[] gradient(double[] values, int from, int to) {
        throw new RuntimeException("Not yet implemented");
    }

    @Override
    public double[] gradientInverse(double[] values, int from, int to) {
        throw new RuntimeException("Not yet implemented");
    }

    @Override
    public String getTransformName() {
        return NAME;
    }

    @Override
    protected double[] transform(double[] values) {
        double propTime = values[0];
        double rateAncestral = Math.exp(values[1]);
        double rateDescendant = Math.exp(values[2]);
        double[] transformed = new double[1];
        transformed[0] = (rateDescendant * propTime) / (rateDescendant * propTime + rateAncestral * (1.0 - propTime));
        return transformed;
    }

    @Override
    protected double[] inverse(double[] values) {
        throw new RuntimeException("Not yet implemented");
    }

    @Override
    protected double getLogJacobian(double[] values) {
        throw new RuntimeException("Not yet implemented");
    }

    @Override
    protected double[] getGradientLogJacobianInverse(double[] values) {
        throw new RuntimeException("Not yet implemented");
    }

    @Override
    public double[][] computeJacobianMatrixInverse(double[] values) {
        throw new RuntimeException("Not yet implemented");
    }

    @Override
    protected boolean isInInteriorDomain(double[] values) {
        // Only the proportion is bounded
        return values[0] >= 0.0 && values[0] <= 1.0;
    }

    public static final AbstractXMLObjectParser PARSER = new AbstractXMLObjectParser() {

        @Override
        public Object parseXMLObject(XMLObject xo) throws XMLParseException {
            final String name = xo.hasId() ? xo.getId() : null;

            TimeToDistanceProportionTransform transform = new TimeToDistanceProportionTransform();

            return transform;
        }

        @Override
        public XMLSyntaxRule[] getSyntaxRules() {
            return new XMLSyntaxRule[0];
        }

        @Override
        public String getParserDescription() {
            return null;
        }

        @Override
        public Class getReturnType() {
            return Transform.MultivariateTransform.class;
        }

        @Override
        public String getParserName() {
            return NAME;
        }
    };
}
