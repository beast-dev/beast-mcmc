/*
 * Copyright (C) 2002-2009 Alexei Drummond and Andrew Rambaut
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
 * BEAST is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with BEAST; if not, write to the
 * Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
 * Boston, MA  02110-1301  USA
 */

package dr.app.beauti.enumTypes;

import dr.evomodelxml.operators.TreeNodeSlideParser;
import dr.inference.operators.ScaleOperator;

/**
 * @author Alexei Drummond
 */
public enum OperatorType {

    SCALE("scale"),
    RANDOM_WALK("randomWalk"),
    RANDOM_WALK_ABSORBING("randomWalkAbsorbing"),
    RANDOM_WALK_REFLECTING("randomWalkReflecting"),
    INTEGER_RANDOM_WALK("integerRandomWalk"),
    UP_DOWN("upDown"),
    UP_DOWN_ALL_RATES_HEIGHTS("upDownAllRatesHeights"),
    SCALE_ALL(ScaleOperator.SCALE_ALL),
    SCALE_INDEPENDENTLY("scaleIndependently"),
    CENTERED_SCALE("centeredScale"),
    DELTA_EXCHANGE("deltaExchange"),
    INTEGER_DELTA_EXCHANGE("integerDeltaExchange"),
    SWAP("swap"),
    BITFLIP("bitFlip"),
    TREE_BIT_MOVE("treeBitMove"),
    SAMPLE_NONACTIVE("sampleNoneActiveOperator"),
    SCALE_WITH_INDICATORS("scaleWithIndicators"),
    UNIFORM("uniform"),
    INTEGER_UNIFORM("integerUniform"),
    SUBTREE_SLIDE("subtreeSlide"),
    NARROW_EXCHANGE("narrowExchange"),
    WIDE_EXCHANGE("wideExchange"),
    GMRF_GIBBS_OPERATOR("gmrfGibbsOperator"),
    WILSON_BALDING("wilsonBalding"),
    NODE_REHIGHT(TreeNodeSlideParser.TREE_NODE_REHEIGHT);

    OperatorType(String displayName) {
        this.displayName = displayName;
    }

    public String toString() {
        return displayName;
    }

    private final String displayName;
}
