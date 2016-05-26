/*
 * TipDateSamplingComponentFactory.java
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

package dr.app.beauti.components.tipdatesampling;

import dr.app.beauti.components.ComponentFactory;
import dr.app.beauti.generator.ComponentGenerator;
import dr.app.beauti.options.ComponentOptions;
import dr.app.beauti.options.BeautiOptions;

/**
 * @author Andrew Rambaut
 * @version $Id$
 */
public class TipDateSamplingComponentFactory implements ComponentFactory {

    private TipDateSamplingComponentFactory() {
        // singleton pattern - private constructor
    }

    @Override
    public Class getOptionsClass() {
        return TipDateSamplingComponentOptions.class;
    }

    public ComponentGenerator createGenerator(final BeautiOptions beautiOptions) {
        return new TipDateSamplingComponentGenerator(beautiOptions);
    }

    public ComponentOptions createOptions(final BeautiOptions beautiOptions) {
        return new TipDateSamplingComponentOptions(beautiOptions);
    }

    public static ComponentFactory INSTANCE = new TipDateSamplingComponentFactory();
}