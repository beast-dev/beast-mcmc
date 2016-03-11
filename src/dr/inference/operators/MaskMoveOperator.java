/*
 * MaskMoveOperator.java
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

import java.util.List;

/**
 * A generic operator that flips masks.
 *
 * @author Marc A. Suchard
 * @version $Id$
 */
public class MaskMoveOperator extends SimpleMCMCOperator {

    public MaskMoveOperator(List<Parameter> masks, Parameter cutPoint, int[] selectBefore, int[] selectAfter, double weight) {
        this.masks = masks;
        this.cutPoint = cutPoint;
        this.selectBefore = selectBefore;
        this.selectAfter = selectAfter;

        setWeight(weight);

        if (!checkMaskValues(masks, cutPoint, selectBefore, selectAfter)) {
            throw new IllegalArgumentException("Bad initialization state");
        }
    }

    public static boolean checkMaskValues(List<Parameter> masks, Parameter cutPoint, int[] selectBefore, int[] selectAfter) {

        int cut = (int) (cutPoint.getParameterValue(0) + 0.5);
        for (int i = 0; i < masks.size(); ++i) {
            Parameter mask = masks.get(i);

            boolean before = (i < cut);

            final int[] ones;
            final int[] zeros;

            if (before) {
                ones = selectBefore;
                zeros = selectAfter;
            } else {
                ones = selectAfter;
                zeros = selectBefore;
            }

            for (int idx : ones) {
                if (mask.getParameterValue(idx) != 1) return false;
            }

            for (int idx: zeros) {
                if (mask.getParameterValue(idx) != 0) return false;
            }
        }
        return true;
    }

    private String printMask() {
        StringBuilder sb = new StringBuilder();
        int cut = (int) (cutPoint.getParameterValue(0) + 0.5);
        sb.append("Cut: " + cut + "\n");
        for (int i = 0; i < masks.size(); ++i) {
            sb.append((i + 1) + " " + masks.get(i).getParameterValue(selectBefore[0]) + " " + masks.get(i).getParameterValue(selectAfter[0]) + "\n");
        }
        return sb.toString();
    }

    public final double doOperation() {
        double logq = 0.0;

        StringBuilder sb = null;

        int currentCutPoint = (int) (cutPoint.getParameterValue(0) + 0.5);

        if (DEBUG) {
            sb = new StringBuilder();
            sb.append("Starting state\n");
            sb.append(printMask());
        }

        final boolean moveUp;
        if (currentCutPoint == 0) {
            moveUp = true;
            logq -= Math.log(2);
        } else if (currentCutPoint == masks.size()) {
            moveUp = false;
            logq -= Math.log(2);
        } else if (MathUtils.nextDouble() < 0.5) {
            moveUp = false;
        } else {
            moveUp = true;
        }

        if (DEBUG) {
            sb.append("moveUp = " + moveUp + "\n");
        }

        Parameter change = (moveUp) ? masks.get(currentCutPoint) : masks.get(currentCutPoint - 1);

        final int[] ones;
        final int[] zeros;
        if (moveUp) {
            ones = selectBefore;
                 zeros = selectAfter;

        } else {


            ones = selectAfter;
                      zeros = selectBefore;
        }

        for (int i : zeros) {
            change.setParameterValueQuietly(i, 0);
        }

        for (int i : ones) {
            change.setParameterValueQuietly(i, 1);
        }

        change.fireParameterChangedEvent();

        cutPoint.setParameterValue(0, currentCutPoint + ((moveUp) ? 1 : -1));

        if (DEBUG) {
            if (!checkMaskValues(masks, cutPoint, selectBefore, selectAfter)) {
                sb.append(printMask());
                System.err.println(sb.toString());
                System.exit(-1);
            }
        }

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
    final private int[] selectBefore;
    final private int[] selectAfter;

    final static private boolean DEBUG = false;
}
