/*
 * JointGradient.java
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

package dr.inference.hmc;

import dr.evomodel.treedatalikelihood.continuous.BranchRateGradient;
import dr.evomodel.treedatalikelihood.discrete.BranchRateGradientForDiscreteTrait;
import dr.xml.*;

import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.*;

/**
 * @author Alexander Fisher
 */

public class JointBranchRateGradient extends JointGradient {

    private static final boolean COMPUTE_IN_PARALLEL = true;

    private final static String JOINT_BRANCH_RATE_GRADIENT = "JointBranchRateGradient";

    public JointBranchRateGradient(List<GradientWrtParameterProvider> derivativeList) {
        super(derivativeList, COMPUTE_IN_PARALLEL ? derivativeList.size() : 0);
    }

    @Override
    double[] getDerivativeLogDensity(DerivativeType derivativeType) {
        return super.getDerivativeLogDensity(derivativeType);
    }

    // **************************************************************
    // XMLObjectParser
    // **************************************************************

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return JOINT_BRANCH_RATE_GRADIENT;
        }

        public Object parseXMLObject(XMLObject xo) {

            List<GradientWrtParameterProvider> derivativeList = new ArrayList<>();

            for (int i = 0; i < xo.getChildCount(); i++) {
                GradientWrtParameterProvider grad = (GradientWrtParameterProvider) xo.getChild(i);
                derivativeList.add(grad);
            }

            return new JointGradient(derivativeList);
        }
        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private final XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
                new OrRule(
                        new ElementRule(BranchRateGradient.class, 1, Integer.MAX_VALUE),
                        new ElementRule(BranchRateGradientForDiscreteTrait.class, 1, Integer.MAX_VALUE)
                ),
        };

        public String getParserDescription() {
            return "Joint branch rate gradient";
        }

        public Class getReturnType() {
            return JointBranchRateGradient.class;
        }
    };
}
