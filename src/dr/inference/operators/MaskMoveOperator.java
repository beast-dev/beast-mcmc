/*
 * MaskFlipOperator.java
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

package dr.inference.operators;

import dr.inference.model.Parameter;
import dr.math.MathUtils;

import java.util.ArrayList;
import java.util.List;

/**
 * A generic operator that flips masks.
 *
 * @author Marc A. Suchard
 * @version $Id$
 */
public class MaskMoveOperator extends SimpleMCMCOperator {

    public MaskMoveOperator(List<Parameter> masks, Parameter cutPoint, Parameter before, Parameter after, double weight) {
        this.masks = masks;
        this.cutPoint = cutPoint;

        selection = new ArrayList<Parameter>(2);
        selection.add(before);
        selection.add(after);

        setWeight(weight);

        if (!checkMaskValues(
                masks, cutPoint, selection.get(0), selection.get(1)
        )) {
            throw new IllegalArgumentException("Bad initialization state");
        }
    }

    public static boolean checkMaskValues(List<Parameter> masks, Parameter cutPoint, Parameter before, Parameter after) {
        for (int i = 0; i < masks.size(); ++i) {
            Parameter mask = masks.get(i);
            Parameter check = (i < (int) (cutPoint.getParameterValue(0) + 0.5)) ? before : after;

            for (int j = 0; j < mask.getDimension(); ++j) {
                if (mask.getParameterValue(j) != check.getParameterValue(j)) {
                    return false;
                }
            }
        }
        return true;
    }

    public final double doOperation() {
        double logq = 0.0;

        int currentCutPoint = (int) (cutPoint.getParameterValue(0) + 0.5);

        final boolean moveUp;
        if (currentCutPoint == 0) {
            moveUp = true;
            logq -= Math.log(2); // TODO Check
        } else if (currentCutPoint == masks.size()) {
            moveUp = false;
            logq -= Math.log(2); // TODO Check
        } else if (MathUtils.nextDouble() < 0.5) {
            moveUp = false;
        } else {
            moveUp = true;
        }

        Parameter change = (moveUp) ? masks.get(currentCutPoint) : masks.get(currentCutPoint - 1);
        Parameter newValues = (moveUp) ? selection.get(1) : selection.get(0);

        for (int i = 0; i < change.getDimension(); ++i) {
            change.setParameterValueQuietly(i, newValues.getParameterValue(i));
        }
        change.fireParameterChangedEvent();

        cutPoint.setParameterValue(0, currentCutPoint + ((moveUp) ? 1 : -1));

        return logq;
    }

    // Interface MCMCOperator
    public final String getOperatorName() {
        return "maskMove(" + cutPoint.getParameterName() + ")";
    }

    public final String getPerformanceSuggestion() {
        return "no performance suggestion";
    }

    public String toString() {
        return getOperatorName();
    }

    // Private instance variables
    final private List<Parameter> masks;
    final private Parameter cutPoint;
    final private List<Parameter> selection;
}
