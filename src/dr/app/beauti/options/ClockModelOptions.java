/*
 * ClockModelOptions.java
 *
 * Copyright (c) 2002-2016 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.tree.UPGMATree;
import dr.evolution.util.Taxa;
import dr.stats.DiscreteStatistics;

import java.util.List;
import java.util.Set;


/**
 * Is this necessary - is likely redundant?
 *
 * @author Alexei Drummond
 * @author Andrew Rambaut
 * @author Walter Xie
 * @version $Id$
 */
public class ClockModelOptions extends ModelOptions {

    private static final long serialVersionUID = 3544930558477534541L;
    // Instance variables
    private final BeautiOptions options;

    public ClockModelOptions(BeautiOptions options) {
        this.options = options;

        initModelParametersAndOpererators();
    }

    /**
     * return a list of parameters that are required
     */

    @Override
    public void initModelParametersAndOpererators() {

    }

    @Override
    public List<Parameter> selectParameters(List<Parameter> params) {
        return null;
    }

    @Override
    public List<Operator> selectOperators(List<Operator> ops) {
        return null;
    }

    @Override
    public String getPrefix() {
        return null;
    }


    public boolean isTipCalibrated() {
        return options.maximumTipHeight > 0;
    }

}
