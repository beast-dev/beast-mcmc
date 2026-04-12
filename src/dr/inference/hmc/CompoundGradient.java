/*
 * CompoundGradient.java
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

package dr.inference.hmc;

import dr.inference.model.*;
import dr.xml.Reportable;

import java.util.ArrayList;
import java.util.IdentityHashMap;
import java.util.List;
import java.util.Map;
import java.util.concurrent.Future;
import java.util.concurrent.atomic.AtomicLong;

/**
 * @author Max Tolkoff
 * @author Marc A. Suchard
 */

public class CompoundGradient implements GradientWrtParameterProvider, HessianWrtParameterProvider, DerivativeWrtParameterProvider, Reportable {

    protected final int dimension;
    final List<GradientWrtParameterProvider> derivativeList;
    private final Likelihood likelihood;
    private final Parameter parameter;

    private final List<DerivativeWrtParameterProvider> newDerivativeList;
    private final DerivativeOrder highestOrder;
    private final ParallelGradientExecutor parallelExecutor;
    private static final AtomicLong DEBUG_SERIAL_CALL_COUNTER = new AtomicLong(0L);

    public CompoundGradient(List<GradientWrtParameterProvider> derivativeList) {
        this(derivativeList, 0);
    }

    public CompoundGradient(List<GradientWrtParameterProvider> derivativeList, int threadCount) {

        this.derivativeList = derivativeList;

        if (derivativeList.size() == 1) {
            likelihood = derivativeList.get(0).getLikelihood();
            parameter = derivativeList.get(0).getParameter();
            dimension = parameter.getDimension();
        } else {
            List<Likelihood> likelihoodList = new ArrayList<>();

            CompoundParameter compoundParameter = new CompoundParameter("hmc") {
                public void fireParameterChangedEvent() {
                    doNotPropagateChangeUp = true;
                    for (Parameter p : uniqueParameters) {
                        p.fireParameterChangedEvent();
                    }
                    doNotPropagateChangeUp = false;
                    fireParameterChangedEvent(-1, ChangeType.ALL_VALUES_CHANGED);
                }
            };

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

            if (Boolean.getBoolean("beast.debug.compoundParameter.identitySummary")) {
                emitCompoundParameterIdentitySummary("CompoundGradient", compoundParameter);
            }
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
            return getDerivativeLogDensitySerialImpl(JointGradient.DerivativeType.GRADIENT);
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

    private double[] getDerivativeLogDensitySerialImpl(JointGradient.DerivativeType derivativeType) {

        double[] result = new double[dimension];
        final boolean debugRepeatability = Boolean.getBoolean("beast.debug.compoundGradientRepeatability");
        final boolean debugContrib = derivativeType == JointGradient.DerivativeType.GRADIENT
                && Boolean.getBoolean("beast.debug.compoundGradientContrib");
        final long debugCall = debugContrib ? DEBUG_SERIAL_CALL_COUNTER.incrementAndGet() : 0L;
        final int debugStride = Integer.getInteger("beast.debug.compoundGradientContribStride", 1);

        int offset = 0;
        int providerIndex = 0;
        for (GradientWrtParameterProvider grad : derivativeList) {

            double[] tmp = derivativeType.getDerivativeLogDensity(grad);
            if (debugRepeatability && derivativeType == JointGradient.DerivativeType.GRADIENT) {
                System.err.println("compoundGradientProvider orderOffset=" + offset
                        + " provider=" + grad.getClass().getSimpleName()
                        + " gradDim=" + grad.getDimension()
                        + " tmpLen=" + tmp.length
                        + " parameter="
                        + (grad.getParameter() == null ? "<null>" : grad.getParameter().getParameterName())
                        + " paramDim="
                        + (grad.getParameter() == null ? -1 : grad.getParameter().getDimension()));
            }
            if (debugRepeatability) {
                final double[] repeat = derivativeType.getDerivativeLogDensity(grad);
                double maxAbs = 0.0;
                int maxIdx = -1;
                for (int i = 0; i < tmp.length; ++i) {
                    final double diff = Math.abs(tmp[i] - repeat[i]);
                    if (diff > maxAbs) {
                        maxAbs = diff;
                        maxIdx = i;
                    }
                }
                if (maxAbs > 0.0) {
                    System.err.println("compoundGradientRepeatability provider=" + grad.getClass().getSimpleName()
                            + " parameter=" + (grad.getParameter() == null ? "<null>" : grad.getParameter().getParameterName())
                            + " maxAbsDiff=" + maxAbs
                            + " maxIdx=" + maxIdx);
                }
            }
            System.arraycopy(tmp, 0, result, offset, grad.getDimension());
            if (debugContrib && debugStride > 0 && (debugCall % debugStride == 0)) {
                emitCompoundGradientDebug(debugCall, providerIndex, grad, offset, tmp);
            }
            offset += grad.getDimension();
            providerIndex++;
        }

        if (debugContrib && debugStride > 0 && (debugCall % debugStride == 0)) {
            emitCompoundGradientDebug(debugCall, -1, null, 0, result);
        }

        return result;
    }

    private static void emitCompoundParameterIdentitySummary(final String label,
                                                             final CompoundParameter compoundParameter) {
        final Map<Parameter, Integer> multiplicity = new IdentityHashMap<Parameter, Integer>();
        for (int i = 0; i < compoundParameter.getParameterCount(); ++i) {
            final Parameter child = compoundParameter.getParameter(i);
            final Integer previous = multiplicity.get(child);
            multiplicity.put(child, previous == null ? 1 : previous + 1);
        }
        int duplicatedChildRefs = 0;
        for (Integer count : multiplicity.values()) {
            if (count != null && count > 1) {
                duplicatedChildRefs++;
            }
        }
        System.err.println("compoundParameterIdentitySummary label=" + label
                + " id=" + System.identityHashCode(compoundParameter)
                + " parameterCount=" + compoundParameter.getParameterCount()
                + " identityUniqueCount=" + multiplicity.size()
                + " duplicatedChildRefs=" + duplicatedChildRefs
                + " dimension=" + compoundParameter.getDimension());
    }

    private void emitCompoundGradientDebug(final long call,
                                           final int providerIndex,
                                           final GradientWrtParameterProvider provider,
                                           final int offset,
                                           final double[] values) {
        final int[] indices = parseDebugIndices();
        final StringBuilder sb = new StringBuilder();
        sb.append("compoundGradientContrib call=").append(call);
        if (providerIndex >= 0) {
            sb.append(" providerIndex=").append(providerIndex);
            sb.append(" providerClass=").append(provider.getClass().getSimpleName());
            sb.append(" parameter=");
            sb.append(provider.getParameter() == null ? "<null>" : provider.getParameter().getParameterName());
            sb.append(" offset=").append(offset);
            for (int globalIdx : indices) {
                final int localIdx = globalIdx - offset;
                if (localIdx >= 0 && localIdx < values.length) {
                    sb.append(" i").append(globalIdx).append("=").append(values[localIdx]);
                } else {
                    sb.append(" i").append(globalIdx).append("=<na>");
                }
            }
        } else {
            sb.append(" providerIndex=sum");
            for (int idx : indices) {
                if (idx >= 0 && idx < values.length) {
                    sb.append(" i").append(idx).append("=").append(values[idx]);
                } else {
                    sb.append(" i").append(idx).append("=<oob>");
                }
            }
        }
        System.err.println(sb.toString());
    }

    private int[] parseDebugIndices() {
        final String raw = System.getProperty("beast.debug.compoundGradient.indices", "6,7,8");
        final String[] tokens = raw.split(",");
        final int[] out = new int[tokens.length];
        for (int i = 0; i < tokens.length; ++i) {
            try {
                out[i] = Integer.parseInt(tokens[i].trim());
            } catch (NumberFormatException nfe) {
                out[i] = -1;
            }
        }
        return out;
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

    @Override
    public double[] getDiagonalHessianLogDensity() {
        if (parallelExecutor != null)
            return getDerivativeLogDensityParallelImpl(JointGradient.DerivativeType.DIAGONAL_HESSIAN);
        else {
            return getDerivativeLogDensitySerialImpl(JointGradient.DerivativeType.DIAGONAL_HESSIAN);
        }
    }

    @Override
    public double[][] getHessianLogDensity() {
        throw new RuntimeException("Not implemented yet");
    }
}
