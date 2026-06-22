/*
 * BMProcessModel.java
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

package dr.evomodel.treedatalikelihood.continuous;

import dr.evomodel.continuous.MultivariateDiffusionModel;
import dr.evomodel.treedatalikelihood.continuous.canonical.math.MatrixOps;
import dr.inference.model.AbstractModel;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;

import java.util.Arrays;

/**
 * Brownian-motion process model with cached trait precision and covariance.
 *
 * <p>The cached snapshot is refreshed only when the wrapped diffusion model
 * changes. Branch-local methods fill caller-owned buffers and scale cached
 * arrays by the effective branch time.
 */
public final class BMProcessModel extends AbstractModel implements ContinuousProcessModel {

    private final int dimension;
    private final MultivariateDiffusionModel diffusionModel;
    private final double[] traitPrecision;
    private final double[] traitCovariance;
    private final double[] cholesky;
    private final double[] inverseWorkspace;

    private double logDetTraitPrecision;
    private boolean dirty = true;

    public BMProcessModel(final MultivariateDiffusionModel diffusionModel) {
        this("bmProcessModel", diffusionModel);
    }

    public BMProcessModel(final String id, final MultivariateDiffusionModel diffusionModel) {
        super(id);
        this.diffusionModel = diffusionModel;
        this.dimension = diffusionModel.getDimension();
        final int dim2 = dimension * dimension;
        this.traitPrecision = new double[dim2];
        this.traitCovariance = new double[dim2];
        this.cholesky = new double[dim2];
        this.inverseWorkspace = new double[dim2];
        addModel(diffusionModel);
        refreshSnapshot();
    }

    @Override
    public int getDimension() {
        return dimension;
    }

    @Override
    public void fillTransitionMoments(final double effectiveTime,
                                      final double[] transitionMatrixOut,
                                      final double[] offsetOut,
                                      final double[] covarianceOut) {
        ensureCurrentSnapshot();
        fillIdentity(transitionMatrixOut);
        Arrays.fill(offsetOut, 0.0);
        for (int i = 0; i < traitCovariance.length; i++) {
            covarianceOut[i] = effectiveTime * traitCovariance[i];
        }
    }

    public void fillTraitPrecisionFlat(final double[] out) {
        ensureCurrentSnapshot();
        System.arraycopy(traitPrecision, 0, out, 0, traitPrecision.length);
    }

    public void fillTraitCovarianceFlat(final double[] out) {
        ensureCurrentSnapshot();
        System.arraycopy(traitCovariance, 0, out, 0, traitCovariance.length);
    }

    public void fillTransitionMatrixFlat(final double[] out) {
        fillIdentity(out);
    }

    public void fillTransitionOffset(final double[] out) {
        Arrays.fill(out, 0.0);
    }

    public void fillTransitionCovarianceFlat(final double effectiveTime, final double[] out) {
        ensureCurrentSnapshot();
        for (int i = 0; i < traitCovariance.length; i++) {
            out[i] = effectiveTime * traitCovariance[i];
        }
    }

    public void fillTransitionPrecision(final double effectiveTime, final double[] out) {
        if (!(effectiveTime > 0.0)) {
            throw new IllegalArgumentException(
                    "BM canonical transitions require a positive effective branch time.");
        }
        ensureCurrentSnapshot();
        final double scale = 1.0 / effectiveTime;
        for (int i = 0; i < traitPrecision.length; i++) {
            out[i] = traitPrecision[i] * scale;
        }
    }

    public double getTransitionLogNormalizer(final double effectiveTime) {
        if (!(effectiveTime > 0.0)) {
            throw new IllegalArgumentException(
                    "BM canonical transitions require a positive effective branch time.");
        }
        ensureCurrentSnapshot();
        return 0.5 * (dimension * MatrixOps.LOG_TWO_PI
                - logDetTraitPrecision
                + dimension * Math.log(effectiveTime));
    }

    public double getTraitPrecisionLogDeterminant() {
        ensureCurrentSnapshot();
        return logDetTraitPrecision;
    }

    public void makeDirty() {
        dirty = true;
    }

    private void refreshSnapshot() {
        MatrixOps.toFlat(diffusionModel.getPrecisionMatrix(), traitPrecision, dimension);
        if (!MatrixOps.tryCholesky(traitPrecision, cholesky, dimension)) {
            throw new IllegalArgumentException("BM diffusion precision matrix is not positive definite.");
        }
        logDetTraitPrecision = MatrixOps.invertFromCholesky(
                cholesky,
                inverseWorkspace,
                traitCovariance,
                dimension);
        dirty = false;
    }

    private void ensureCurrentSnapshot() {
        if (dirty) {
            refreshSnapshot();
        }
    }

    private void fillIdentity(final double[] out) {
        Arrays.fill(out, 0.0);
        for (int i = 0; i < dimension; i++) {
            out[i * dimension + i] = 1.0;
        }
    }

    @Override
    protected void handleModelChangedEvent(final Model model, final Object object, final int index) {
        dirty = true;
    }

    @Override
    protected void handleVariableChangedEvent(final Variable variable,
                                              final int index,
                                              final Parameter.ChangeType type) {
        dirty = true;
    }

    @Override
    protected void storeState() {
    }

    @Override
    protected void restoreState() {
        dirty = true;
    }

    @Override
    protected void acceptState() {
    }
}
