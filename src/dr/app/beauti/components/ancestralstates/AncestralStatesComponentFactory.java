/*
 * AncestralStatesComponentFactory.java
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

package dr.app.beauti.components.ancestralstates;

import dr.app.beauti.components.ComponentFactory;
import dr.app.beauti.generator.ComponentGenerator;
import dr.app.beauti.options.BeautiOptions;
import dr.app.beauti.options.ComponentOptions;

/**
 * @author Andrew Rambaut
 */
public class AncestralStatesComponentFactory implements ComponentFactory {

    private AncestralStatesComponentFactory() {
        // singleton pattern - private constructor
    }

    @Override
    public Class getOptionsClass() {
        return AncestralStatesComponentOptions.class;
    }

    @Override
    public ComponentGenerator createGenerator(final BeautiOptions beautiOptions) {
        return new AncestralStatesComponentGenerator(beautiOptions);
    }

    @Override
    public ComponentOptions createOptions(final BeautiOptions beautiOptions) {
        return new AncestralStatesComponentOptions();
    }

    public static ComponentFactory INSTANCE = new AncestralStatesComponentFactory();
}
