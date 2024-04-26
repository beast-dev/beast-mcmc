/*
 * DemographicModel.java
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

package dr.evomodel.coalescent.demographicmodel;

import dr.evolution.coalescent.DemographicFunction;
import dr.evolution.util.Units;
import dr.inference.model.AbstractModel;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;

/**
 * This interface provides methods that describe a demographic model.
 * <p/>
 * Parts of this class were derived from C++ code provided by Oliver Pybus.
 *
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @author Korbinian Strimmer
 * @version $Id: DemographicModel.java,v 1.28 2005/09/26 14:27:38 rambaut Exp $
 */
public abstract class DemographicModel extends AbstractModel implements Units {

    /**
     * abstract base class for parametric demographic models.
     * @param name name of the XML element
     */
    public DemographicModel(String name) {
        this(name, 0.0);
    }

    /**
     * abstract base class for parametric demographic models.
     * @param name name of the XML element
     * @param timeOffset an offset in time scale for use when multiple demographic models are being used
     */
    public DemographicModel(String name, double timeOffset) {
        super(name);

        this.timeOffset = timeOffset;
    }

    // general functions

    public abstract DemographicFunction getDemographicFunction();

    // **************************************************************
    // Model IMPLEMENTATION
    // **************************************************************

    protected void handleModelChangedEvent(Model model, Object object, int index) {
        // no intermediates need to be recalculated...
    }

    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        // no intermediates need to be recalculated...
    }

    protected void storeState() {
    } // no additional state needs storing

    protected void restoreState() {
    } // no additional state needs restoring

    protected void acceptState() {
    } // no additional state needs accepting

    public double getTimeOffset() {
        return timeOffset;
    }

    public void setTimeOffset(double timeOffset) {
        this.timeOffset = timeOffset;
    }

    private double timeOffset;

    // **************************************************************
    // Units IMPLEMENTATION
    // **************************************************************

    /**
     * Units in which population size is measured.
     */
    private Type units;

    /**
     * sets units of measurement.
     *
     * @param u units
     */
    public void setUnits(Type u) {
        units = u;
    }

    /**
     * returns units of measurement.
     */
    public Type getUnits() {
        return units;
    }
}