/*
 * CompoundGradient.java
 *
 * Copyright (c) 2002-2017 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

import dr.inference.model.*;
import dr.xml.Reportable;

import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.Future;

/**
 * @author Max Tolkoff
 * @author Marc A. Suchard
 */

public class CompoundGradient implements GradientWrtParameterProvider, DerivativeWrtParameterProvider, Reportable {

    protected final int dimension;
    final List<GradientWrtParameterProvider> derivativeList;
    private final Likelihood likelihood;
    private final Parameter parameter;

    private final List<DerivativeWrtParameterProvider> newDerivativeList;
    private final DerivativeOrder highestOrder;
    private final ParallelGradientExecutor parallelExecutor;

    CompoundGradient(List<GradientWrtParameterProvider> derivativeList) {
        this(derivativeList, 0);
    }

    CompoundGradient(List<GradientWrtParameterProvider> derivativeList, int threadCount) {

        this.derivativeList = derivativeList;

        if (derivativeList.size() == 1) {
            likelihood = derivativeList.get(0).getLikelihood();
            parameter = derivativeList.get(0).getParameter();
            dimension = parameter.getDimension();
        } else {
            List<Likelihood> likelihoodList = new ArrayList<>();
            CompoundParameter compoundParameter = new CompoundParameter("hmc");

            int dim = 0;
            for (GradientWrtParameterProvider grad : derivativeList) {
                for (Likelihood likelihood : grad.getLikelihood().getLikelihoodSet()) {
                    if (!(likelihoodList.contains(likelihood))) {
                        likelihoodList.add(likelihood);
                    }
                }

                Parameter p = grad.getParameter();
                compoundParameter.addParameter(p);

                dim += p.getDimension();
            }

            likelihood = new CompoundLikelihood(likelihoodList);
            parameter = compoundParameter;
            dimension = dim;
        }

        // NEW
        this.newDerivativeList = new ArrayList<>();
        for (GradientWrtParameterProvider p : derivativeList) {
            if (p instanceof DerivativeWrtParameterProvider) { // TODO Remove if
                DerivativeWrtParameterProvider provider = (DerivativeWrtParameterProvider) p;
                newDerivativeList.add(provider);
            }
        }
        this.highestOrder = DerivativeWrtParameterProvider.getHighestOrder(newDerivativeList);

        // Parallel threading

        if (threadCount > 1 || threadCount < 0) {
            parallelExecutor = new ParallelGradientExecutor(threadCount, derivativeList);
        } else {
            parallelExecutor = null;
        }

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
    public int getDimension(DerivativeOrder order) {
        return order.getDerivativeDimension(dimension);
    }

    @Override
    public int getDimension() {
        return dimension;
    }

    @Override
    public double[] getDerivativeLogDensity(DerivativeOrder order) {

        assert (highestOrder.getValue() >= order.getValue());

        double[] result = new double[dimension];

          int offset = 0;
          for (DerivativeWrtParameterProvider provider : newDerivativeList) {

              double[] tmp = provider.getDerivativeLogDensity(order);
              System.arraycopy(tmp, 0, result, offset, tmp.length);
              offset += tmp.length;
          }

          return result;
    }

    @Override
    public DerivativeOrder getHighestOrder() {
        return highestOrder;
    }
    
    @Override
    public double[] getGradientLogDensity() {
        if (parallelExecutor != null)
            return getDerivativeLogDensityParallelImpl(JointGradient.DerivativeType.GRADIENT);
        else {
            return getDerivativeLogDensitySerialImpl();
        }
    }

    private double[] getDerivativeLogDensityParallelImpl(JointGradient.DerivativeType derivativeType) {

        return parallelExecutor.getDerivativeLogDensityInParallel(derivativeType, (gradients, length) -> {
            double[] reduction = new double[length];
            int offset = 0;
            for (Future<double[]> result : gradients) {
                double[] tmp = result.get();
                System.arraycopy(tmp, 0, reduction, offset, tmp.length);
                offset += tmp.length;
            }
            return reduction;
        }, dimension);
    }

    private double[] getDerivativeLogDensitySerialImpl() {

        double[] result = new double[dimension];

        int offset = 0;
        for (GradientWrtParameterProvider grad : derivativeList) {
            
            double[] tmp = grad.getGradientLogDensity(); // TODO Generalize for Hessian
            System.arraycopy(tmp, 0, result, offset, grad.getDimension());
            offset += grad.getDimension();
        }

        return result;
    }

    @Override
    public String getReport() {
        return  "compoundGradient." + parameter.getParameterName() + "\n" +
                GradientWrtParameterProvider.getReportAndCheckForError(this,
                Double.NEGATIVE_INFINITY, Double.POSITIVE_INFINITY,
                GradientWrtParameterProvider.TOLERANCE);
    }

    public List<GradientWrtParameterProvider> getDerivativeList() {
        return derivativeList;
    }
}
