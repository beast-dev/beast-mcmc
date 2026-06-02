/*
 * TransformedParameterOperator.java
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

package dr.inference.operators;

import dr.inference.model.BoundedSpace;
import dr.inference.model.Parameter;
import dr.inference.model.TransformedParameter;

public class TransformedParameterOperator extends AbstractAdaptableOperator {
    private final boolean isAdaptable;
    private final SimpleMCMCOperator subOperator;
    private final TransformedParameter parameter;
    private final boolean checkValid;
    private final BoundedSpace generalBounds;
    public static final String TRANSFORMED_OPERATOR = "transformedParameterOperator";

    public TransformedParameterOperator(SimpleMCMCOperator operator, BoundedSpace generalBounds) {

        this.subOperator = operator;
        setWeight(operator.getWeight());
        this.isAdaptable = operator instanceof AbstractAdaptableOperator;
        this.parameter = (TransformedParameter) operator.getParameter();

        this.generalBounds = generalBounds;
        this.checkValid = generalBounds != null;
    }


    @Override
    protected void setAdaptableParameterValue(double value) {
        if (isAdaptable) {
            ((AbstractAdaptableOperator) subOperator).setAdaptableParameterValue(value);
        }
    }

    @Override
    protected double getAdaptableParameterValue() {
        if (isAdaptable) {
            return ((AbstractAdaptableOperator) subOperator).getAdaptableParameterValue();
        }
        return 0;
    }

    @Override
    public double getRawParameter() {
        if (isAdaptable) {
            return ((AbstractAdaptableOperator) subOperator).getRawParameter();
        }
        throw new RuntimeException("not actually adaptable parameter");
    }

    @Override
    public String getAdaptableParameterName() {
        if (isAdaptable) {
            return ((AbstractAdaptableOperator) subOperator).getAdaptableParameterName();
        }
        throw new RuntimeException("not actually adaptable parameter");
    }

    @Override
    public String getOperatorName() {
        return TRANSFORMED_OPERATOR + "." + subOperator.getOperatorName();
    }

    @Override
    public double doOperation() {
        double[] oldValues = parameter.getParameterUntransformedValues();
        double ratio = subOperator.doOperation();
        double[] newValues = parameter.getParameterUntransformedValues();


        if (checkValid) { // GH: below is sloppy, but best I could do without refactoring how Parameter handles bounds
            if (generalBounds == null && !parameter.isWithinBounds()) {
                return Double.NEGATIVE_INFINITY;
            } else if (!generalBounds.isWithinBounds(parameter.getParameterValues())) {
                return Double.NEGATIVE_INFINITY;
            }
        }

        // Compute Jacobians
        ratio += parameter.diffLogJacobian(oldValues, newValues);

        return ratio;
    }

    @Override
    public Parameter getParameter() {
        return subOperator.getParameter();
    }
}
