/*
 * OperatorUtils.java
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

package dr.inference.operators;

/**
 * Utility functions to aid in operator optimization etc.
 *
 * @author Alexei Drummond
 * @version $Id: OperatorUtils.java,v 1.3 2005/05/24 20:26:00 rambaut Exp $
 */
public class OperatorUtils {

    /**
     * Given a current window size and an acceptance level and target acceptance level,
     * returns a new window size.
     *
     * @param delta        the current delta (half the window size)
     * @param currentLevel the current acceptance probability
     * @param targetLevel  the target acceptance probability
     * @return a new window size
     */
    public static double optimizeWindowSize(double delta, double currentLevel, double targetLevel) {
        return optimizeWindowSize(delta, Double.MAX_VALUE, currentLevel, targetLevel);
    }

    /**
     * Given a current delta (+-) and an acceptance level and target acceptance level,
     * returns a new delta.
     *
     * @param delta        the current delta (half the window size)
     * @param maxDelta     the maximum delta allowed
     * @param currentLevel the current acceptance probability
     * @param targetLevel  the target acceptance probability
     * @return a new window size
     */
    public static double optimizeWindowSize(double delta, double maxDelta, double currentLevel, double targetLevel) {

        if (delta <= 0.0) {
            throw new IllegalArgumentException("random walk window size cannot be negative: " + delta);
        }

        double ratio = currentLevel / targetLevel;

        if (ratio > 2.0) ratio = 2.0;
        if (ratio < 0.5) ratio = 0.5;

        double newDelta = delta * ratio;


        if (newDelta > maxDelta) newDelta = maxDelta;

        return newDelta;
    }

    /**
     * Given a current scale factor and acceptance level and a target acceptance level,
     * returns a new scale factor.
     *
     * @param scaleFactor  the current scale factor
     * @param currentLevel the current acceptance probability
     * @param targetLevel  the target acceptance probability
     * @return a new scale factor
     */
    public static double optimizeScaleFactor(double scaleFactor, double currentLevel, double targetLevel) {

        if (scaleFactor <= 0.0 || scaleFactor >= 1.0) {
            throw new IllegalArgumentException("scale factor was " + scaleFactor + "!");
        }

        double ratio = currentLevel / targetLevel;

        if (ratio > 2.0) ratio = 2.0;
        if (ratio < 0.5) ratio = 0.5;

        // new scale factor
        return Math.pow(scaleFactor, ratio);
	}
}
