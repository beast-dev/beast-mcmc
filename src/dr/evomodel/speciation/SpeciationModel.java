/*
 * SpeciationModel.java
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

package dr.evomodel.speciation;

import dr.evolution.tree.Tree;
import dr.evolution.util.Taxon;
import dr.evolution.util.Units;
import dr.inference.model.AbstractModel;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;

import java.util.Set;

/**
 * @author Alexei Drummond
 * @author Andrew Rambaut
 */
public abstract class SpeciationModel extends AbstractModel implements Units {
    /**
     * Units in which population size is measured.
     */
    private Units.Type units;

    public SpeciationModel(String name, Type units) {
        super(name);
        setUnits(units);
    }

    public abstract double calculateTreeLogLikelihood(Tree tree);

    public abstract double calculateTreeLogLikelihood(Tree tree, Set<Taxon> exclude);

    // True if Yule.
    //
    // Not abstract - non supporting derived classes do not need to override anything
    public boolean isYule() {
        return false;
    }

    // Likelihood for the speciation model conditional on monophyly and calibration densities in
    // 'calibration'.
    //
    // The likelihood enforces the monophyly, so there is no need to specify it again in the XML.
    //
    public double calculateTreeLogLikelihood(Tree tree, CalibrationPoints calibration) {
        return Double.NEGATIVE_INFINITY;
    }

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

    /**
     * sets units of measurement.
     *
     * @param u units
     */
    public void setUnits(Units.Type u) {
        units = u;
    }

    /**
     * returns units of measurement.
     */
    public Units.Type getUnits() {
        return units;
    }
}
