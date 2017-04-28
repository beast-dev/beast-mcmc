/*
 * AppendedPotentialDerivativeParser.java
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

import dr.inference.distribution.DistributionLikelihood;
import dr.inference.distribution.MultivariateDistributionLikelihood;
import dr.inference.hmc.CompoundGradient;
import dr.inference.hmc.GradientWrtParameterProvider;
import dr.inference.model.*;
import dr.xml.*;

import java.util.ArrayList;
import java.util.List;

/**
 * @author Max Tolkoff
 */
public class AppendedPotentialDerivativeParser extends AbstractXMLObjectParser {

    public final static String SUM_DERIVATIVE = "appendedPotentialDerivative";
    public static final String SUM_DERIVATIVE2 = "compoundGradient";

    @Override
    public String getParserName() {
        return SUM_DERIVATIVE;
    }

    @Override
    public String[] getParserNames() {
        return new String[] { SUM_DERIVATIVE, SUM_DERIVATIVE2 };
    }

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        List<GradientWrtParameterProvider> gradList = new ArrayList<GradientWrtParameterProvider>();
        List<Likelihood> likelihoodList = new ArrayList<Likelihood>(); // TODO Remove?

        for (int i = 0; i < xo.getChildCount(); ++i) {

            Object obj = xo.getChild(i);

            GradientWrtParameterProvider grad;
            Likelihood likelihood;

            if (obj instanceof DistributionLikelihood) {
                DistributionLikelihood dl = (DistributionLikelihood) obj;
                if (!(dl.getDistribution() instanceof GradientProvider)) {
                    throw new XMLParseException("Not a gradient provider");
                }

                throw new RuntimeException("Not yet implemented");

            } else if (obj instanceof MultivariateDistributionLikelihood) {
                final MultivariateDistributionLikelihood mdl = (MultivariateDistributionLikelihood) obj;
                if (!(mdl.getDistribution() instanceof GradientProvider)) {
                    throw new XMLParseException("Not a gradient provider");
                }

                final GradientProvider provider = (GradientProvider) mdl.getDistribution();
                final Parameter parameter = mdl.getDataParameter();
                likelihood = mdl;

                grad = new GradientWrtParameterProvider.ParameterWrapper(provider, parameter, mdl);

            } else if (obj instanceof GradientWrtParameterProvider) {
                grad = (GradientWrtParameterProvider) obj;
                likelihood = grad.getLikelihood();
            } else {
                throw new XMLParseException("Not a Gaussian process");
            }

            gradList.add(grad);
            likelihoodList.add(likelihood);
        }

        return new CompoundGradient(gradList);
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(GradientWrtParameterProvider.class, 1, Integer.MAX_VALUE),
    };

    @Override
    public String getParserDescription() {
        return null;
    }

    @Override
    public Class getReturnType() {
        return CompoundGradient.class;
    }
}
