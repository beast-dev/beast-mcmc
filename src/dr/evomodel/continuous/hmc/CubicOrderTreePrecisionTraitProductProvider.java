/*
 * CubicOrderTreePrecisionTraitProductProvider.java
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

package dr.evomodel.continuous.hmc;

import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.evomodel.treedatalikelihood.continuous.ContinuousDataLikelihoodDelegate;
import dr.inference.model.Model;
import dr.inference.model.Parameter;

/**
 * @author Marc A. Suchard
 */
public class CubicOrderTreePrecisionTraitProductProvider extends TreePrecisionTraitProductProvider {

    public CubicOrderTreePrecisionTraitProductProvider(TreeDataLikelihood treeDataLikelihood,
                                                       ContinuousDataLikelihoodDelegate likelihoodDelegate) {
        super(treeDataLikelihood, likelihoodDelegate);
        this.isPrecisionKnown = false;
    }

    @Override
    public double[] getProduct(Parameter vector) {
        return expensiveProduct(vector, getTreeTraitPrecision());
    }

    @Override
    public double[] getMassVector() {

        double[][] treeTraitVariance = likelihoodDelegate.getTreeTraitVariance();
        final int dim = treeTraitVariance.length;

        double[] mass = new double[dim];
        for (int i = 0; i < dim; ++i) {
            mass[i] = treeTraitVariance[i][i];
        }

        return mass;
    }

    @Override
    public double getTimeScale() {


        double[][] treeTraitVariance = likelihoodDelegate.getTreeTraitVariance();
        final int dim = treeTraitVariance.length;

        double max = Double.MIN_VALUE;
        for (int i = 0; i < dim; ++i) {
            max = Math.max(max, treeTraitVariance[i][i]);
        }

        return Math.sqrt(max);
    }

    @Override
    public double getTimeScaleEigen() {
        return 0;
    }

    private double[][] getTreeTraitPrecision() {

        if (!DO_CACHE) {
            treeTraitPrecision = likelihoodDelegate.getTreeTraitPrecision();
        } else if (!isPrecisionKnown) {
            treeTraitPrecision = likelihoodDelegate.getTreeTraitPrecision();
            isPrecisionKnown = true;
        }

        return treeTraitPrecision;
    }

    @Override
    protected void handleModelChangedEvent(Model model, Object object, int index) {
        if (model == tree) {
            isPrecisionKnown = false;
        }
    }

    @Override
    protected void storeState() {
        if (DO_CACHE) {
            final int dim = treeTraitPrecision.length;
            if (savedTreeTraitPrecision == null) {
                savedTreeTraitPrecision = new double[dim][dim];
            }

            for (int i = 0; i < dim; ++i) {
                System.arraycopy(treeTraitPrecision[i], 0, savedTreeTraitPrecision[i], 0, dim);
            }

            savedIsPrecisionKnown = isPrecisionKnown;
        }
    }

    @Override
    protected void restoreState() {
        if (DO_CACHE) {
            double[][] tmp = treeTraitPrecision;
            treeTraitPrecision = savedTreeTraitPrecision;
            savedTreeTraitPrecision = tmp;

            isPrecisionKnown = savedIsPrecisionKnown;
        }
    }

    private static final boolean DO_CACHE = true;

    private boolean isPrecisionKnown;
    private boolean savedIsPrecisionKnown;
    private double[][] treeTraitPrecision;
    private double[][] savedTreeTraitPrecision;
}
