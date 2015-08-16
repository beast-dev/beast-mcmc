/*
 * MigrationModel.java
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

package dr.evomodel.coalescent.structure;

import dr.evolution.colouring.ColourChangeMatrix;
import dr.inference.model.AbstractModel;
import dr.inference.model.Model;

/**
 * @version $Id: MigrationModel.java,v 1.3 2006/09/08 14:28:07 rambaut Exp $
 *
 * @author Andrew Rambaut
 * @author Alexei Drummond
 */
public abstract class MigrationModel extends AbstractModel
{

	public MigrationModel(String name) { super(name); }

	// general functions

	public abstract ColourChangeMatrix getMigrationMatrix();

    public abstract double[] getMigrationRates(double time);

    // **************************************************************
    // Model IMPLEMENTATION
    // **************************************************************

	protected void handleModelChangedEvent(Model model, Object object, int index) {
		// no intermediates need to be recalculated...
	}

	protected void storeState() {} // no additional state needs storing
	protected void restoreState() {} // no additional state needs restoring
	protected void acceptState() {} // no additional state needs accepting

}