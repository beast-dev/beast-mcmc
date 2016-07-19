/*
 * AbstractPopulationMRCAModel.java
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

package dr.oldevomodel.approxPopTree;

import dr.evolution.alignment.Patterns;
import dr.evolution.tree.NodeRef;
import dr.inference.model.AbstractModel;

import java.util.LinkedList;

/**
 * Package: AbstractPopulationMRCAModel
 * Description:
 * <p/>
 * <p/>
 * Created by
 *
 * @author Alexander V. Alekseyenko (alexander.alekseyenko@gmail.com)
 *         Date: Feb 1, 2010
 *         Time: 5:30:01 PM
 */
abstract public class AbstractPopulationMRCAModel extends AbstractModel {
    protected double time; // time of the population node above
    protected double tMRCA;

    public AbstractPopulationMRCAModel(String name, double populationTime) {
        super(name);
        time = populationTime;
    }

    public double getMRCATime(LinkedList<NodeRef> nodes) { // get time of actual MRCA for the sequence nodes below
        return tMRCA;
    }

    abstract public double drawMRCATime(LinkedList<NodeRef> nodes); // re-draw from the distribution of t(MRCA)

    abstract public double[][] getMRCAPartials(LinkedList<NodeRef> nodes, Patterns patterns); // get partials at Population node

    abstract public double[][] drawMRCAPartials(LinkedList<NodeRef> nodes, Patterns patterns); // re-draw partial at Population node

    protected double[][] computeProfilePartials(LinkedList<NodeRef> nodes, Patterns patterns) { // compute profile partials for the sequence nodes, useful for the above methods

        int patternCount = patterns.getPatternCount();
        int stateCount = patterns.getStateCount();
        double partials[][] = new double[patternCount][stateCount];
        // Do the profile computation here
        //
        return partials;
    }
}
