/*
 * LoadingsHamiltonianMCParser.java
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

import dr.inference.operators.hmc.deprecated.LoadingsHamiltonianMC;
import dr.inference.distribution.MomentDistributionModel;
import dr.inference.model.LatentFactorModel;
import dr.inference.model.MatrixParameter;
import dr.inference.operators.CoercionMode;
import dr.xml.*;

/**
 * Created by max on 1/11/16.
 */
@Deprecated
public class LoadingsHamiltonianMCParser extends AbstractXMLObjectParser {
    public static final String LOADINGS_HAMILTONIAN_MC="loadingsHamiltonianMC";
    public static final String WEIGHT="weight";
    public static final String STEP_SIZE="stepSize";
    public static final String N_STEPS="nSteps";
    public static final String MOMENTUM_SD="momentumSd";

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        LatentFactorModel lfm=(LatentFactorModel) xo.getChild(LatentFactorModel.class);
        MomentDistributionModel prior=(MomentDistributionModel) xo.getChild(MomentDistributionModel.class);
        double weight=xo.getDoubleAttribute(WEIGHT);
        CoercionMode mode=CoercionMode.parseMode(xo);
        int nSteps=xo.getIntegerAttribute(N_STEPS);
        double stepSize=xo.getDoubleAttribute(STEP_SIZE);
        double momentumSd= xo.getDoubleAttribute(MOMENTUM_SD);
        MatrixParameter loadings=(MatrixParameter) xo.getChild(MatrixParameter.class);

        return new LoadingsHamiltonianMC(lfm, prior, weight, mode, stepSize, nSteps, momentumSd, loadings);
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
            new ElementRule(MomentDistributionModel.class),
            new ElementRule(MatrixParameter.class),
    };


    @Override
    public String getParserDescription() {
        return "Hamiltonian Monte Carlo for loadings matrix in a latent factor model";
    }

    @Override
    public Class getReturnType() {
        return LoadingsHamiltonianMC.class;
    }

    @Override
    public String getParserName() {
        return LOADINGS_HAMILTONIAN_MC;
    }
}
