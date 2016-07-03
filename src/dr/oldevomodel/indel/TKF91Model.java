/*
 * TKF91Model.java
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

package dr.oldevomodel.indel;

import dr.oldevomodelxml.indel.TKF91ModelParser;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;
import org.w3c.dom.Document;
import org.w3c.dom.Element;

/**
 * This class models a constant mutation rate
 * (parameter: mu = mutation rate). <BR>
 *
 * @author Alexei Drummond
 * @version $Id: TKF91Model.java,v 1.7 2005/05/24 20:25:57 rambaut Exp $
 */
public class TKF91Model extends IndelModel {

    private final Parameter lengthDistParameter, deathRateParameter;

    public TKF91Model(Parameter lengthDistParameter, Parameter deathRateParameter, Type units) {
        super(TKF91ModelParser.TKF91_MODEL);

        this.lengthDistParameter = lengthDistParameter;
        addVariable(lengthDistParameter);

        this.deathRateParameter = deathRateParameter;
        addVariable(deathRateParameter);

        setUnits(units);
    }

    public final double getLengthDistributionValue() {
        return lengthDistParameter.getParameterValue(0);
    }

    public final double getBirthRate(int length) {

        throw new RuntimeException("Not implemented");
        //return birthRateParameter.getParameterValue(0);

    }

    public final double getDeathRate(int length) {
        return deathRateParameter.getParameterValue(0);
    }

    // *****************************************************************
    // Interface ModelComponent
    // *****************************************************************

    public String getModelComponentName() {
        return TKF91ModelParser.TKF91_MODEL;
    }

    protected void handleModelChangedEvent(Model model, Object object, int index) {
        // Substitution model has changed so fire model changed event
        listenerHelper.fireModelChanged(this, object, index);
    }

    protected final void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        // no intermediates need to be recalculated...
    }

    protected void storeState() {
    } // no extra state apart from parameters

    protected void acceptState() {
    } // no extra state apart from parameters

    protected void restoreState() {
    } // no extra state apart from parameters

    // **************************************************************
    // XMLElement IMPLEMENTATION
    // **************************************************************

    public Element createElement(Document doc) {
        throw new RuntimeException("Not implemented!");
    }

}
