/*
 * MaskableSpeciationModel.java
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

package dr.evomodel.speciation;

/**
 * Created by IntelliJ IDEA.
 * User: adru001
 * Date: Nov 2, 2010
 * Time: 10:56:03 AM
 * To change this template use File | Settings | File Templates.
 */
public abstract class MaskableSpeciationModel extends SpeciationModel {

    public MaskableSpeciationModel(String name, Type units) {
        super(name, units);
    }

    // a model specific implementation that allows this speciation model
    // to be partially masked by another -- useful in model averaging applications
    public abstract void mask(SpeciationModel mask);

    public abstract void unmask();
}
