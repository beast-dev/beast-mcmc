/*
 * SpeciationModelGradientProvider.java
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

package dr.evomodel.speciation;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.inference.model.Parameter;

public interface SpeciationModelGradientProvider {

    default double getNodeHeightGradient(Tree tree, NodeRef node) {
        throw new RuntimeException("Not yet implemented");
    }

    default double[] getBirthRateGradient(Tree tree, NodeRef node) {
        throw new RuntimeException("Not yet implemented");
    }

    default double[] getDeathRateGradient(Tree tree, NodeRef node) {
        throw new RuntimeException("Not yet implemented");
    }

    default double[] getSamplingRateGradient(Tree tree, NodeRef node) {
        throw new RuntimeException("Not yet implemented");
    }

    default double[] getTreatmentProbabilityGradient(Tree tree, NodeRef node) {
        throw new RuntimeException("Not yet implemented");
    }

    default double[] getSamplingProbabilityGradient(Tree tree, NodeRef node) {
        throw new RuntimeException("Not yet implemented");
    }

    default Parameter getBirthRateParameter() {
        throw new RuntimeException("Not yet implemented");
    }

    default Parameter getDeathRateParameter() {
        throw new RuntimeException("Not yet implemented");
    }

    default Parameter getSamplingRateParameter() {
        throw new RuntimeException("Not yet implemented");
    }

    default Parameter getTreatmentProbabilityParameter() {
        throw new RuntimeException("Not yet implemented");
    }

    default Parameter getSamplingProbabilityParameter() {
        throw new RuntimeException("Not yet implemented");
    }

    default int getBirthRateIndex() { throw new RuntimeException("Not implemented"); }

    default int getDeathRateIndex() { throw new RuntimeException("Not implemented"); }

    default int getSamplingRateIndex() { throw new RuntimeException("Not implemented"); }

    default int getSamplingProbabilityIndex() { throw new RuntimeException("Not implemented"); }

    default int getTreatmentProbabilityIndex() { throw new RuntimeException("Not implemented"); }

    default double[] getBreakPoints() { throw new RuntimeException("Not yet implemented"); }
    
    default void precomputeGradientConstants() { throw new RuntimeException("Not yet implemented"); }

    default void processGradientModelSegmentBreakPoint(double[] gradient,
                                                       int currentModelSegment,
                                                       double intervalStart,
                                                       double segmentIntervalEnd,
                                                       int nLineages) {
        throw new RuntimeException("Not yet implemented");
    }

    default void processGradientInterval(double[] gradient,
                                         int currentModelSegment,
                                         double intervalStart,
                                         double intervalEnd,
                                         int nLineages) { throw new RuntimeException("Not yet implemented"); }

    default void processGradientSampling(double[] gradient,
                                         int currentModelSegment,
                                         double intervalEnd) { throw new RuntimeException("Not yet implemented"); }

    default void processGradientCoalescence(double[] gradient,
                                            int currentModelSegment,
                                            double intervalEnd) { throw new RuntimeException("Not yet implemented"); }

    default void processGradientOrigin(double[] gradient,
                                       int currentModelSegment,
                                       double totalDuration) { throw new RuntimeException("Not yet implemented"); }

    default void logConditioningProbability(int currentModelSegment, double[] gradient) { throw new RuntimeException("Not yet implemented"); }

    default void updateGradientModelValues(int currentModelSegment) { throw new RuntimeException("Not yet implemented"); }

    default int getGradientLength() { throw new RuntimeException("Not yet implemented"); }
}
