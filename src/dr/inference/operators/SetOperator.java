/*
 * SetOperator.java
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

package dr.inference.operators;

import dr.inference.model.Parameter;
import dr.math.MathUtils;
import org.w3c.dom.Document;
import org.w3c.dom.Element;

/**
 * A generic operator for selecting uniformly from a discrete set of values.
 *
 * @author Alexei Drummond
 * @author Andrew Rambaut
 * @version $Id: SetOperator.java,v 1.12 2005/05/24 20:26:00 rambaut Exp $
 */
public class SetOperator extends SimpleMCMCOperator {

    public SetOperator(Parameter parameter, double[] values) {
        this.parameter = parameter;
        this.values = values;
    }

    /**
     * @return the parameter this operator acts on.
     */
    public Parameter getParameter() {
        return parameter;
    }

    /**
     * change the parameter and return the hastings ratio.
     */
    public final double doOperation() {

        int index = MathUtils.nextInt(values.length);
        double newValue = values[index];

        if (newValue < parameter.getBounds().getLowerLimit(index) || newValue > parameter.getBounds().getUpperLimit(index)) {
//            throw new OperatorFailedException("proposed value outside boundaries");
            return Double.NEGATIVE_INFINITY;
        }

        parameter.setParameterValue(index, newValue);

        return 0.0;
    }

    public Element createOperatorElement(Document document) {
        throw new RuntimeException("Not implememented!");
    }

    public String getOperatorName() {
        return "setOperator(" + parameter.getParameterName() + ")";
    }

    public String getPerformanceSuggestion() {
        return "No suggestions";
    }

    //PRIVATE STUFF

    private Parameter parameter = null;
    private double[] values;

}
