/*
 * EqualRatesProvider.java
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

package dr.evomodel.substmodel;

import dr.evomodel.tree.TreeModel;

/**
 * Rate provider for the equal-rates model on any data type: Q_ij = pi_j / (1 - sum_i pi_i^2) (JC69, F81, the Mk
 * model, the binary model), or Q_ij = pi_j unnormalised. Every pair has the same exchangeability, so the model
 * is hierarchical on every character tree with lambda_root = 1 / (1 - sum pi^2) (or 1) and lambda_m = 0 at every
 * other internal node. lambda_root > 0 always, so the generator is always valid and there is nothing to bound.
 * <p>
 * The provider can also take an existing GeneralSubstitutionModel (binarySubstitutionModel,
 * generalSubstitutionModel, including the BSSVS form) whose relative rates are all equal; it reads the model's
 * generator after every change and refuses unequal rates.
 *
 * @author Xiang Ji
 */
public class EqualRatesProvider extends AbstractHierarchicalRateProvider {

    public static final double DEFAULT_TOLERANCE = 1e-10;

    private final boolean normalised;
    private final GeneralSubstitutionModel source;
    private final double[] q;

    public EqualRatesProvider(TreeModel tree, FrequencyModel frequency, boolean normalised) {
        this("equalRates", tree, frequency, normalised, null);
    }

    private EqualRatesProvider(String name, TreeModel tree, FrequencyModel frequency, boolean normalised,
                               GeneralSubstitutionModel source) {
        super(name, tree, frequency);
        this.normalised = normalised;
        this.source = source;
        this.q = source == null ? null : new double[stateCount * stateCount];
        if (source != null) {
            addModel(source);
            getNormalization();
        }
    }

    public static EqualRatesProvider of(GeneralSubstitutionModel model, TreeModel tree) {
        if (model instanceof ComplexSubstitutionModel) {
            throw new IllegalArgumentException("ratesFrom needs a reversible model; "
                    + model.getClass().getSimpleName() + " is not");
        }
        return new EqualRatesProvider(model.getId() == null ? "equalRates" : model.getId() + ".equalRates", tree,
                model.getFrequencyModel(), true, model);
    }

    @Override
    protected double computeLambdas(double[] lambda) {
        final int root = tree.getRoot().getNumber();
        if (source != null) {
            source.getInfinitesimalMatrix(q);
            double s = Double.NaN;
            for (int i = 0; i < stateCount; i++) {
                for (int j = 0; j < stateCount; j++) {
                    if (i != j) {
                        final double value = q[i * stateCount + j] / frequency.getFrequency(j);
                        if (Double.isNaN(s)) {
                            s = value;
                        } else if (Math.abs(value - s) > DEFAULT_TOLERANCE * Math.max(1.0, Math.abs(s))) {
                            throw new IllegalArgumentException(source.getModelName() + " does not have equal rates ("
                                    + dataType.getCode(i) + "-" + dataType.getCode(j) + ": " + value + " vs " + s
                                    + "), so it is not the equal-rates model");
                        }
                    }
                }
            }
            lambda[root] = s;
            return 1.0;
        }
        double sumSquares = 0.0;
        for (int i = 0; i < stateCount; i++) {
            final double f = frequency.getFrequency(i);
            sumSquares += f * f;
        }
        final double norm = 1.0 - sumSquares;
        lambda[root] = normalised ? 1.0 / norm : 1.0;
        return normalised ? norm : 1.0;
    }
}
