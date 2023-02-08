/*
 * HessianWrapperParser.java
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
import dr.inference.hmc.HessianWrtParameterProvider;
import dr.inference.model.HessianProvider;
import dr.inference.model.Parameter;
import dr.xml.*;

/**
 * @author Marc A. Suchard
 * @author Xiang Ji
 */

public class HessianWrapperParser  extends AbstractXMLObjectParser {

    public static final String NAME = "hessian";
    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        Object obj = xo.getChild(0);

        if (obj instanceof MultivariateDistributionLikelihood) {
            final MultivariateDistributionLikelihood mdl = (MultivariateDistributionLikelihood)
                    xo.getChild(MultivariateDistributionLikelihood.class);

            if (!(mdl.getDistribution() instanceof HessianProvider)) {
                throw new XMLParseException("Not a hessian provider");
            }

            final HessianProvider provider = (HessianProvider) mdl.getDistribution();
            final Parameter parameter = mdl.getDataParameter();

            return new HessianWrtParameterProvider.ParameterWrapper(provider, parameter, mdl);
        } else if (obj instanceof DistributionLikelihood) {
            DistributionLikelihood dl = (DistributionLikelihood) obj;
            if (!(dl.getDistribution() instanceof HessianProvider)) {
                throw new XMLParseException("Not a hessian provider");
            }

            final HessianProvider provider = (HessianProvider) dl.getDistribution();
            final Parameter parameter = (Parameter) xo.getChild(Parameter.class);

            // TODO Ensure that parameter and data inside provider are the same.

            return new HessianWrtParameterProvider.ParameterWrapper(provider, parameter, dl);
        } else {
            assert(obj instanceof HessianWrtParameterProvider);
            return obj;
        }
    }


    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return new XMLSyntaxRule[] {
                new XORRule(
                        new XORRule(
                                new ElementRule(MultivariateDistributionLikelihood.class),
                                new AndRule(
                                        new ElementRule(DistributionLikelihood.class),
                                        new ElementRule(Parameter.class)
                                )
                        ),
                        new ElementRule(HessianWrtParameterProvider.class)
                )
        };
    }

    @Override
    public String getParserDescription() {
        return null;
    }

    @Override
    public Class getReturnType() {
        return HessianWrtParameterProvider.class;
    }

    @Override
    public String getParserName() {
        return NAME;
    }
}
