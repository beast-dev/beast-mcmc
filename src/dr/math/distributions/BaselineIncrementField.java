/*
 * BaselineIncrementField.java
 *
 * Copyright (c) 2002-2023 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.math.distributions;

import dr.inference.distribution.RandomField;
import dr.inference.model.GradientProvider;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;

/**
 * @author Marc Suchard
 * @author Yucai Shao
 * @author Andy Magee
 */
public class BaselineIncrementField extends RandomFieldDistribution {

    public static final String TYPE = "BaselineIncrementField";

    private final Parameter baseline;
    private final Parameter increments;
    private final RandomField.WeightProvider weights;

    public BaselineIncrementField(String name,
                                  Parameter baseline,
                                  Parameter increments,
                                  RandomField.WeightProvider weights) {
        super(name);

        this.baseline = baseline;
        this.increments = increments;
        this.weights = weights;

        addVariable(baseline);
        addVariable(increments);

        if (weights != null) {
            addModel(weights);
            throw new IllegalArgumentException("Unsure how weights influence this field");
        }
    }

    @Override
    public Variable<Double> getLocationVariable() {
        throw new RuntimeException("Not yet implemented");
    }

    @Override
    public double[] nextRandom() {
        throw new RuntimeException("Not yet implemented");
    }

    @Override
    protected void handleModelChangedEvent(Model model, Object object, int index) {

    }

    @Override
    protected void storeState() {

    }

    @Override
    protected void restoreState() {

    }

    @Override
    protected void acceptState() {

    }

    @Override
    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {

    }

    @Override
    public int getDimension() {
        throw new RuntimeException("Not yet implemented");
    }

    @Override
    public double[] getGradientLogDensity(Object x) {
        throw new RuntimeException("Not yet implemented");
    }

    @Override
    public double[] getDiagonalHessianLogDensity(Object x) {
        throw new RuntimeException("Not yet implemented");
    }

    @Override
    public double[][] getHessianLogDensity(Object x) {
        throw new RuntimeException("Not yet implemented");
    }

    @Override
    public double logPdf(double[] x) {
        throw new RuntimeException("Not yet implemented");
    }

    @Override
    public double[][] getScaleMatrix() {
        throw new RuntimeException("Not yet implemented");
    }

    @Override
    public double[] getMean() {
        throw new RuntimeException("Not yet implemented");
    }

    @Override
    public String getType() { return TYPE; }

    @Override
    public GradientProvider getGradientWrt(Parameter parameter) {
        throw new RuntimeException("Not yet implemented");
    }
}
