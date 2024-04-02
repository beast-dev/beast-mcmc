/*
 * ContinuousTraitPartialsProvider.java
 *
 * Copyright (c) 2002-2019 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.evomodel.treedatalikelihood.continuous;

import dr.evolution.tree.Tree;
import dr.evomodel.treedatalikelihood.continuous.cdi.PrecisionType;
import dr.evomodel.treedatalikelihood.preorder.WrappedNormalSufficientStatistics;
import dr.inference.model.CompoundParameter;

import java.util.ArrayList;
import java.util.List;

/**
 * @author Marc A. Suchard
 */

public interface ContinuousTraitPartialsProvider {

    boolean bufferTips();

    int getTraitCount();

    int getTraitDimension();

    String getTipTraitName();

    void setTipTraitName(String name);

    default int getDataDimension() {
        return getTraitDimension();
    }

    PrecisionType getPrecisionType();

    double[] getTipPartial(int taxonIndex, boolean fullyObserved);

    @Deprecated
    List<Integer> getMissingIndices(); // use getTraitMissingIndicators() instead

    boolean[] getDataMissingIndicators(); // returns null for no missing data

    default boolean[] getTraitMissingIndicators() { // returns null for no missing traits
        return getDataMissingIndicators();
    }

    CompoundParameter getParameter();

    String getModelName();

    boolean usesMissingIndices();

    ContinuousTraitPartialsProvider[] getChildModels();

    default double[] drawTraitsBelowConditionalOnDataAndTraitsAbove(double[] aboveTraits) {
        throw new RuntimeException("Conditional sampling not yet implemented for " + this.getClass());
    }

    default double[] transformTreeTraits(double[] traits) {
        return traits;
    }

    default boolean getDefaultAllowSingular() {
        return false;
    }

    default boolean suppliesWishartStatistics() {
        return true;
    }

    default int[] getPartitionDimensions() {
        return new int[]{getTraitDimension()};
    }

    default void addTreeAndRateModel(Tree treeModel, ContinuousRateTransformation rateTransformation) {
        // Do nothing
    }

    default WrappedNormalSufficientStatistics partitionNormalStatistics(WrappedNormalSufficientStatistics statistic,
                                                                        ContinuousTraitPartialsProvider provider) {
        if (this == provider) {
            return statistic;
        }
        throw new RuntimeException("This class does not currently support 'partitionNormalStatistics' with " +
                "a provider other than itself.");
    }

    default ContinuousTraitPartialsProvider getProviderForTrait(String trait) {
        if (trait.equals(getTipTraitName())) {
            return this;
        }
        throw new RuntimeException("Partials provider does not have trait '" + trait + "'");
    }

    static boolean[] indicesToIndicator(List<Integer> indices, int n) {

        if (indices == null) {
            return null;
        }

        boolean[] indicator = new boolean[n];

        for (int i : indices) {
            indicator[i] = true;
        }

        return indicator;

    }

    static List<Integer> indicatorToIndices(boolean[] indicators) { //TODO: test
        List<Integer> indices = new ArrayList<>();

        for (int i = 0; i < indicators.length; i++) {
            if (indicators[i]) {
                indices.add(i);
            }
        }

        return indices;
    }

}
