/*
 * SiteModelOptions.java
 *
 * Copyright (c) 2002-2015 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.app.beauti.options;

import dr.app.beauti.types.OperatorType;

import java.util.List;

/**
 * @author Walter Xie
 */
public class SiteModelOptions extends ModelOptions {
    private static final long serialVersionUID = -3347506415688390314L;

    // Instance variables
    private final BeautiOptions options;


    public SiteModelOptions(BeautiOptions options) {
        this.options = options;
    }


    @Override
    public void initModelParametersAndOpererators() {
        createOperator("dataLikelihoodMultivariate", "Multiple", "Adaptive Multivariate Normal", "multiple",
                OperatorType.ADAPTIVE_MULTIVARIATE, 1.0, treeWeights);
        createOperator("treePriorMultivariate", "Multiple", "Adaptive Multivariate Normal", "multiple",
                OperatorType.ADAPTIVE_MULTIVARIATE, 1.0, treeWeights);

    }

    @Override
    public List<Parameter> selectParameters(List<Parameter> params) {
        return params;
    }

    @Override
    public List<Operator> selectOperators(List<Operator> ops) {
        ops.add(getOperator("dataLikelihoodMultivariate"));
        ops.add(getOperator("treePriorMultivariate"));
        return ops;
    }

    @Override
    public String getPrefix() {
        return null;
    }

}
