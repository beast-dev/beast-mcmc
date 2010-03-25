/*
 * Copyright (C) 2002-2009 Alexei Drummond and Andrew Rambaut
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
 * BEAST is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with BEAST; if not, write to the
 * Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
 * Boston, MA  02110-1301  USA
 */

package dr.app.beauti.options;


import dr.app.beauti.enumTypes.OperatorType;
import dr.app.beauti.enumTypes.PriorScaleType;
import dr.evomodel.substmodel.AbstractSubstitutionModel;
import dr.evomodel.substmodel.SVSGeneralSubstitutionModel;
import dr.inference.operators.RateBitExchangeOperator;

import java.util.List;

/**
 * @author Walter Xie
 */
public class DiscreteTraitOptions extends TraitsOptions {

    public static enum LocationSubstModelType {
        SYM_SUBST("Symmetric substitution model"),
        ASYM_SUBST("Asymmetric substitution model");

        LocationSubstModelType(String name) {
            this.name = name;
        }

        public String toString() {
            return name;
        }

        private final String name;
    }
    
    public static final String PHYLOGEOGRAPHIC = "phylogeographic ";

    private LocationSubstModelType locationSubstType = LocationSubstModelType.SYM_SUBST;
    private boolean activateBSSVS = false;

    public DiscreteTraitOptions(TraitGuesser traitGuesser) {
        super(traitGuesser);
    }

    @Override
    protected void initTraitParametersAndOperators() {

        createParameterUniformPrior("frequencies", PHYLOGEOGRAPHIC + "base frequencies", PriorScaleType.UNITY_SCALE, 0.25, 0.0, 1.0);
        createCachedGammaPrior("rates", "location substitution model rates",
                PriorScaleType.SUBSTITUTION_PARAMETER_SCALE, 1.0, 1.0, 1.0, 0, Double.POSITIVE_INFINITY, false);
        createParameter("indicators", "location substitution model rate indicators");

        createParameterExponentialPrior("mu", PHYLOGEOGRAPHIC + "mutation rate parameter",
                PriorScaleType.SUBSTITUTION_PARAMETER_SCALE, 0.1, 1.0, 0.0, 0.0, 10.0);

        createDiscreteStatistic("nonZeroRates", "for mutation rate parameter (if BSSVS was selected)");  // BSSVS was selected

        createOperator("rates", OperatorType.SCALE_INDEPENDENTLY, demoTuning, 30);
        createScaleOperator("mu", "for mutation rate parameter", demoTuning, 10);

        createOperator("indicators", OperatorType.BITFLIP, -1.0, 30);
        createTagInsideOperator(OperatorType.BITFIP_IN_SUBST.toString(), "mu",
                "bit Flip In Substitution Model Operator", "mu", OperatorType.BITFIP_IN_SUBST,
                SVSGeneralSubstitutionModel.SVS_GENERAL_SUBSTITUTION_MODEL, PREFIX_ + AbstractSubstitutionModel.MODEL, demoTuning, 30);
        // <svsGeneralSubstitutionModel idref="originModel"/>
        createOperatorUsing2Parameters(RateBitExchangeOperator.OPERATOR_NAME, "(indicators, rates)",
                "rateBitExchangeOperator (If both BSSVS and asymmetric subst selected)",
                "indicators", "rates", OperatorType.RATE_BIT_EXCHANGE, -1.0, 6.0);
    }


    /**
     * return a list of parameters that are required
     *
     * @param params the parameter list
     */
    @Override
    public void selectParameters(List<Parameter> params) {
        params.add(getParameter("frequencies"));
        params.add(getParameter("rates"));
        params.add(getParameter("mu"));

        if (activateBSSVS) {
            params.add(getParameter("indicators"));         
            params.add(getParameter("nonZeroRates"));
        }
    }

    /**
     * return a list of operators that are required
     *
     * @param ops the operator list
     */
    @Override
    public void selectOperators(List<Operator> ops) {
        ops.add(getOperator("rates"));
        ops.add(getOperator("mu"));

        if (activateBSSVS) {
            ops.add(getOperator("indicators"));
            ops.add(getOperator(OperatorType.BITFIP_IN_SUBST.toString()));

            if (locationSubstType == LocationSubstModelType.ASYM_SUBST)
                ops.add(getOperator(RateBitExchangeOperator.OPERATOR_NAME));
        }

    }

    @Override
    public boolean isSpecifiedTraitAnalysis(String traitName) {
        return this.partitionName == traitName;
    }

    /////////////////////////////////////////////////////////////

    public LocationSubstModelType getLocationSubstType() {
        return locationSubstType;
    }

    public void setLocationSubstType(LocationSubstModelType locationSubstType) {
        this.locationSubstType = locationSubstType;
    }

    public boolean isActivateBSSVS() {
        return activateBSSVS;
    }

    public void setActivateBSSVS(boolean activateBSSVS) {
        this.activateBSSVS = activateBSSVS;
    }

}

