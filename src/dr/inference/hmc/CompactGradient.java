/*
 * CompactGradient.java
 *
 * Copyright (c) 2002-2022 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.inference.hmc;

import dr.inference.model.CompoundParameter;
import dr.inference.model.Likelihood;
import dr.inference.model.Parameter;
import dr.xml.Reportable;

/**
 * @author Marc A. Suchard
 */
public class CompactGradient implements HessianWrtParameterProvider, Reportable {

    private final GradientWrtParameterProvider source;
    private final Parameter sourceParameter;

    private final Likelihood likelihood;
    private final Parameter parameter;
    private final int[] map;
    private final int dimension;

    public CompactGradient(GradientWrtParameterProvider source) {

        this.source = source;
        this.sourceParameter = source.getParameter();
        this.likelihood = source.getLikelihood();

        ParameterMap map = constructParameter(sourceParameter);
        this.parameter = map.parameter;
        this.map = map.map;
        this.dimension = parameter.getDimension();
    }

    @Override
    public Likelihood getLikelihood() {
        return likelihood;
    }

    @Override
    public Parameter getParameter() {
        return parameter;
    }

    @Override
    public int getDimension() {
        return dimension;
    }

    @Override
    public double[] getGradientLogDensity() {

        double[] input = source.getGradientLogDensity();
        return compact(input);
    }

    @Override
    public double[] getDiagonalHessianLogDensity() {

        if (!(source instanceof HessianWrtParameterProvider)) {
            throw new RuntimeException("Must use Hessian providers");
        }

        double[] input = ((HessianWrtParameterProvider) source).getDiagonalHessianLogDensity();
        return compact(input);
    }

    @Override
    public double[][] getHessianLogDensity() {
        throw new RuntimeException("Not yet implemented");
    }

    @Override
    public String getReport() {
        return  "compactGradient." + sourceParameter.getParameterName() + "\n" +
                GradientWrtParameterProvider.getReportAndCheckForError(this,
                        Double.NEGATIVE_INFINITY, Double.POSITIVE_INFINITY,
                        GradientWrtParameterProvider.TOLERANCE);
    }

    private double[] compact(double[] input) {
        double[] output = new double[dimension];

        for (int i = 0; i < map.length; ++i) {
            output[map[i]] += input[i];
        }

        return output;
    }

    private void map(Parameter p, int offset0, int offset1, CompoundParameter destination, int[] map) {
        if (destination != null) {
            destination.addParameter(p);
        }

        for (int k = 0; k < p.getDimension(); ++k) {
            int from = offset0 + k;
            int to = offset1 + k;
            map[from] = to;
        }
    }

    private static class ParameterMap {
        CompoundParameter parameter;
        int[] map;

        ParameterMap(int dim) {
            parameter = new CompoundParameter("compact");
            map = new int[dim];
        }
    }

    private ParameterMap constructParameter(Parameter source) {

        ParameterMap map = new ParameterMap(sourceParameter.getDimension());

        if (source instanceof CompoundParameter) {
            CompoundParameter cp = (CompoundParameter) source;

            int current = 0;
            for (int i = 0; i < cp.getParameterCount(); ++i) {
                Parameter p = cp.getParameter(i);

                int past = 0;
                boolean found = false;
                for (int j = 0; j < i && !found; ++j) {
                    Parameter q = cp.getParameter(j);
                    if (p == q) {
                        map(p, current, past, null, map.map);
                        found = true;
                    }
                    past += q.getDimension();
                }
                if (!found) {
                    map(p, current, current, map.parameter, map.map);
                }

                current += p.getDimension();
            }
        } else {
            throw new IllegalArgumentException("Can only compact compound gradients");
        }

        return map;
    }
}
