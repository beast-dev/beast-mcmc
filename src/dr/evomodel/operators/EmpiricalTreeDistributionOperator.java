/*
 * EmpiricalTreeDistributionOperator.java
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

package dr.evomodel.operators;

import dr.inference.operators.SimpleMCMCOperator;
import dr.evomodel.tree.EmpiricalTreeDistributionModel;

/**
 * @author Andrew Rambaut
 * @version $Id$
 */
public class EmpiricalTreeDistributionOperator extends SimpleMCMCOperator {

    public static final String EMPIRICAL_TREE_DISTRIBUTION_OPERATOR = "empiricalTreeDistributionOperator";
    public static final String METROPOLIS_HASTINGS = "metropolisHastings";

    private final EmpiricalTreeDistributionModel treeModel;
    private final boolean metropolisHastings;

    public EmpiricalTreeDistributionOperator(EmpiricalTreeDistributionModel treeModel, boolean metropolisHastings, double weight) {
        this.treeModel = treeModel;
        setWeight(weight);
        this.metropolisHastings = metropolisHastings;
    }

    // IMPLEMENTATION: SimpleMCMCOperator
    // =========================================

    public double doOperation() {
        treeModel.drawTreeIndex();

        if (metropolisHastings) {           
            return 0.0;
        }

        // if not MetropolisHastings, always accept these moves...
        return Double.POSITIVE_INFINITY;
    }

    public String getPerformanceSuggestion() {
        return "";
    }

    public String getOperatorName() {
        return EMPIRICAL_TREE_DISTRIBUTION_OPERATOR;
    }

}