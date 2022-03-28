/*
 * BranchRates.java
 *
 * Copyright (c) 2002-2020 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.evomodel.branchratemodel;

import dr.evolution.tree.BranchRates;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.inference.model.Parameter;

import java.util.function.DoubleBinaryOperator;

/**
 * @author Marc A. Suchard
 * @author Alexander Fisher
 */
public interface DifferentiableBranchRates  extends BranchRates {

    default Tree getTree() { return null; } // TODO Deprecate

    default double getUntransformedBranchRate(Tree tree, NodeRef node) { return getBranchRate(tree, node); }

    double getBranchRateDifferential(Tree tree, NodeRef node); // TODO Deprecate

    double getBranchRateSecondDifferential(Tree tree, NodeRef node); // TODO Deprecate

    Parameter getRateParameter();

    int getParameterIndexFromNode(NodeRef node); // TODO deprecate

    ArbitraryBranchRates.BranchRateTransform getTransform(); // TODO Should remove

    double[] updateGradientLogDensity(double[] gradient, double[] value,
                                      int from, int to);

    double[] updateDiagonalHessianLogDensity(double[] diagonalHessian, double[] gradient, double[] value,
                                             int from, int to);

    default double mapReduceOverRates(NodeRateMap map, DoubleBinaryOperator reduce, double initial) {
        throw new RuntimeException("Not implemented"); // TODO Suggests this should be in subclass
    }

    default void forEachOverRates(NodeRateMap map) {
        throw new RuntimeException("Not implemented"); // TODO Suggests this should be in subclass
    }
}
