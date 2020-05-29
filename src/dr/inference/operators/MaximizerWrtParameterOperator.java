/*
 * MaximizerWrtParameterOperator.java
 *
 * Copyright (c) 2002-2016 Alexei Drummond, Andrew Rambaut and Marc Suchard
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


package dr.inference.operators;

import dr.evomodel.treedatalikelihood.discrete.MaximizerWrtParameter;
import dr.inferencexml.operators.MaximizerWrtParameterOperatorParser;

/**
 * @author Marc Suchard
 * @author Xiang Ji
 */
public class MaximizerWrtParameterOperator extends SimpleMCMCOperator implements GibbsOperator {

    private final MaximizerWrtParameter maximizerWrtParameter;

    private final int maximumOperations;

    private int numOperations = 0;

    public MaximizerWrtParameterOperator(MaximizerWrtParameter maximizerWrtParameter,
                                         int maximumOperations,
                                         double weight) {
        this.maximizerWrtParameter = maximizerWrtParameter;
        this.maximumOperations = maximumOperations;
        setWeight(weight);
    }

    @Override
    public String getOperatorName() {
        return MaximizerWrtParameterOperatorParser.MAXIMIZER_PARAMETER_OPERATOR;
    }

    @Override
    public double doOperation() {
        if (numOperations < maximumOperations) {
            maximizerWrtParameter.maximize();
        }
        numOperations++;
        return 0;
    }
}
