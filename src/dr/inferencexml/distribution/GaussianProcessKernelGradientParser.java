/*
 * GaussianProcessKernelParser.java
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

package dr.inferencexml.distribution;

import dr.inference.distribution.RandomField;
import dr.inference.model.Parameter;
import dr.math.distributions.gp.AdditiveGaussianProcessDistribution;
import dr.math.distributions.gp.BasisDimension;
import dr.math.distributions.gp.GaussianProcessKernel;
import dr.math.distributions.gp.GaussianProcessKernelGradient;
import dr.xml.*;

import java.util.List;

public class GaussianProcessKernelGradientParser extends AbstractXMLObjectParser {

    private static final String PARSER_NAME = "gaussianProcessKernelGradient";

    public String getParserName() { return PARSER_NAME; }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        RandomField randomField = (RandomField) xo.getChild(RandomField.class);
        List<Parameter> parameters;

        if (randomField.getDistribution() instanceof AdditiveGaussianProcessDistribution) {
            BasisDimension basis = null;
            AdditiveGaussianProcessDistribution distribution = (AdditiveGaussianProcessDistribution) randomField.getDistribution();
            if (xo.getChild(GaussianProcessKernel.class) != null) {
                GaussianProcessKernel kernel = (GaussianProcessKernel) xo.getChild(GaussianProcessKernel.class);
                if (kernel.getParameters().isEmpty()) {
                    throw new XMLParseException("The kernel has no parameters");
                }
                parameters = kernel.getParameters();
                for (BasisDimension b : distribution.getBases()) {
                    if (b.getKernel() == kernel) {
                        basis = b;
                        break;
                    }
                }
                if (basis == null) {
                    throw new XMLParseException("The kernel is not consistent with the random field");
                }
            } else {
                parameters = xo.getAllChildren(Parameter.class);
                if (parameters.isEmpty()) {
                    throw new XMLParseException("No parameters are specified");
                } else {
                    for (BasisDimension b : distribution.getBases()) {
                        if (b.getKernel().getParameters().containsAll(parameters)) {
                            basis = b;
                            List<Parameter> kernelParameters = basis.getKernel().getParameters();
                            if (!kernelParameters.containsAll(parameters)) {
                                throw new XMLParseException("The chosen parameters are inconsistent with the kernel");
                            }
                            break;
                        }
                    }
                    if (basis == null) {
                        throw new XMLParseException("The parameters must belong to the same kernel " +
                                "within one basis of the distribution in the random field");
                    }
                }
            }
            return new GaussianProcessKernelGradient(randomField, basis, parameters);
        } else {
            throw new XMLParseException("Not yet implemented for this distribution");
        }
    }

    public XMLSyntaxRule[] getSyntaxRules() { return rules; }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(RandomField.class),
            new XORRule(
                    new XMLSyntaxRule[]{
                        new ElementRule(GaussianProcessKernel.class),
                        new ElementRule(Parameter.class, 1, 2),
                    }
            ),
    };

    public String getParserDescription() {
        return PARSER_NAME;
    }

    public Class getReturnType() { return GaussianProcessKernelGradient.class; }
}
