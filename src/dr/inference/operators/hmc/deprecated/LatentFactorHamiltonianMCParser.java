/*
 * LatentFactorHamiltonianMCParser.java
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

package dr.inference.operators.hmc.deprecated;

import dr.evomodel.continuous.FullyConjugateMultivariateTraitLikelihood;
import dr.inference.operators.hmc.deprecated.LatentFactorHamiltonianMC;
import dr.inference.model.LatentFactorModel;
import dr.inference.operators.CoercionMode;
import dr.xml.*;

/**
 * Created by max on 12/2/15.
 */
@Deprecated
public class LatentFactorHamiltonianMCParser extends AbstractXMLObjectParser {
    public static final String LATENT_FACTOR_MODEL_HAMILTONIAN_MC="LatentFactorHamiltonianMC";
    public static final String WEIGHT="weight";
    public static final String N_STEPS="nSteps";
    public static final String STEP_SIZE="stepSize";
    public static final String MOMENTUM_SD="momentumSd";
    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        LatentFactorModel lfm=(LatentFactorModel) xo.getChild(LatentFactorModel.class);
        FullyConjugateMultivariateTraitLikelihood tree=(FullyConjugateMultivariateTraitLikelihood) xo.getChild(FullyConjugateMultivariateTraitLikelihood.class);
        double weight=xo.getDoubleAttribute(WEIGHT);
        CoercionMode mode=CoercionMode.parseMode(xo);
        int nSteps=xo.getIntegerAttribute(N_STEPS);
        double stepSize=xo.getDoubleAttribute(STEP_SIZE);
        double momentumSd= xo.getDoubleAttribute(MOMENTUM_SD);


        return new LatentFactorHamiltonianMC(lfm, tree, weight, mode, stepSize, nSteps, momentumSd);


    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private static final XMLSyntaxRule[] rules = {
            AttributeRule.newDoubleRule(WEIGHT),
            AttributeRule.newDoubleRule(STEP_SIZE),
            AttributeRule.newIntegerRule(N_STEPS),
            AttributeRule.newDoubleRule(MOMENTUM_SD),
            new ElementRule(LatentFactorModel.class),
            new ElementRule(FullyConjugateMultivariateTraitLikelihood.class),
    };

    @Override
    public String getParserDescription() {
        return "Hamiltonian Monte Carlo For factors";
    }

    @Override
    public Class getReturnType() {
        return LatentFactorHamiltonianMC.class;
    }

    @Override
    public String getParserName() {
        return LATENT_FACTOR_MODEL_HAMILTONIAN_MC;
    }
}
