/*
 * ExperimentalGeneralizedLinearModelParser.java
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

package dr.inferencexml.glm;

import dr.inference.distribution.GeneralizedLinearModel;
import dr.inference.distribution.LogisticRegression;
import dr.inference.glm.ExperimentalGeneralizedLinearModel;
import dr.inference.glm.LogisticRegressionDesignMatrixGradient;
import dr.inference.glm.LogisticRegressionGradientWrtParameter;
import dr.inference.glm.LogisticRegressionRandomEffectsGradientWrtParameter;
import dr.inference.hmc.GradientWrtParameterProvider;
import dr.inference.model.Parameter;
import dr.math.distributions.SplineBasisMatrix;
import dr.xml.*;

import static dr.inferencexml.operators.hmc.HamiltonianMonteCarloOperatorParser.MASK;

/**
 * @author Marc A Suchard
 * @author Filippo Monti
 */
public class LogisticRegressionGradientParser extends AbstractXMLObjectParser {

    public static final String GLM_GRADIENT = "glmGradient";

    public String getParserName() {
        return GLM_GRADIENT;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        GeneralizedLinearModel glm = (GeneralizedLinearModel) xo.getChild(GeneralizedLinearModel.class);
        Parameter parameter = (Parameter) xo.getChild(Parameter.class);

        final GradientWrtParameterProvider gradient;
        int index = glm.getEffectNumber(parameter);
        if (index >= 0) {
            gradient = new LogisticRegressionGradientWrtParameter(
                    (LogisticRegression) glm, parameter);
        } else {
            index = -1;
            for (int i = 0; i < glm.getNumberOfRandomEffects(); i++) {
                if (glm.getRandomEffect(i) == parameter) {
                    index = i;
                    break;
                }
            }
            if (index >= 0) {
                gradient = new LogisticRegressionRandomEffectsGradientWrtParameter(
                        (LogisticRegression) glm, parameter);
            } else {
                for (int i = 0; i < glm.getNumberOfFixedEffects(); ++i) {
                    if (glm.getDesignMatrix(i) == parameter) {
                        index = i;
                        break;
                    }
                }
                if (index >= 0) {
                    gradient = new LogisticRegressionDesignMatrixGradient(
                            (LogisticRegression) glm, (SplineBasisMatrix) parameter);
                } else {
                    throw new XMLParseException("Unable to find independent parameter '" + parameter.getId() +
                            "' in '" + xo.getId() + "'");
                }
            }
        }

        return gradient;
    }

    public XMLSyntaxRule[] getSyntaxRules() { return rules; }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(GeneralizedLinearModel.class),
            new ElementRule(Parameter.class),
            new ElementRule(MASK, new XMLSyntaxRule[] {
                    new ElementRule(Parameter.class),
            }, true),
    };

    public String getParserDescription() {
        return "Calculates the generalized linear model likelihood of the dependent parameters given one or more blocks of independent parameters and their design matrix.";
    }

    public Class getReturnType() {
        return ExperimentalGeneralizedLinearModel.class;
    }
}
