/*
 * CompoundGaussianProcess.java
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

package dr.math.distributions;

import dr.inference.distribution.DistributionLikelihood;
import dr.inference.model.CompoundLikelihood;
import dr.inference.model.Likelihood;
import dr.xml.Reportable;

import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.*;

/**
 * @author Marc A. Suchard
 */

public class CompoundGaussianProcess implements GaussianProcessRandomGenerator, Reportable {

    private final List<GaussianProcessRandomGenerator> gpList;
    private final List<Integer> copyList;
    private final List<Likelihood> likelihoodList;
    private final CompoundLikelihood compoundLikelihood;

    private final ExecutorService pool;
    private final int threadCount;
    private final List<Callable<DrawResult>> callers;

    private static final boolean USE_POOL = false;

    private int dimension = -1;

    public CompoundGaussianProcess(List<GaussianProcessRandomGenerator> gpList, List<Likelihood> likelihoodList,
                                   List<Integer> copyList) {
        this.gpList = gpList;
        this.copyList = copyList;
        this.likelihoodList = likelihoodList;
        compoundLikelihood = new CompoundLikelihood(likelihoodList);

        if (USE_POOL) {
            callers = createTasks();
            threadCount = callers.size();
            pool = Executors.newFixedThreadPool(threadCount);
        } else {
            callers = null;
            threadCount = -1;
            pool = null;
        }
    }

    public boolean contains(Likelihood likelihood) {
        return likelihoodList.contains(likelihood);
    }

    public int getDimension() {
        if (dimension == -1) {
            dimension = 0;
            for (GaussianProcessRandomGenerator gp : gpList) {
                dimension += gp.getDimension();
            }
        }
        return dimension;
    }

    @Override
    public double[][] getPrecisionMatrix() {
        if (gpList.size() == 1) {
            return gpList.get(0).getPrecisionMatrix();
        } else {
            final int dim = getDimension();
            double[][] precision = new double[dim][dim];

            int offset = 0;
            for (GaussianProcessRandomGenerator gp : gpList) {
                final int d = gp.getDimension();
                double[][] p = gp.getPrecisionMatrix();

                for (int i = 0; i < d; ++i) {
                    System.arraycopy(p[i], 0, precision[offset + i], offset, d);
                }

                offset += d;
            }

            return precision;
        }
    }

    @Override
    public Likelihood getLikelihood() { return compoundLikelihood; }

    @Override
    public String getReport() {
        StringBuilder sb = new StringBuilder();
        sb.append("compoundGP: " + getLikelihood().getLogLikelihood());
        return sb.toString();
    }

    private class DrawResult {
        final double[] result;
        final int offset;

        DrawResult(double[] result, int offset) {
            this.result = result;
            this.offset = offset;
        }
    }

    private class DrawCaller implements Callable<DrawResult> {

        public DrawCaller(GaussianProcessRandomGenerator gp, int copies, int offset, boolean isUnivariate) {
            this.gp = gp;
            this.copies = copies;
            this.offset = offset;
            this.isUnivariate = isUnivariate;
        }

        public DrawResult call() throws Exception {

            final double[] vector;
            if (isUnivariate) {
                vector = new double[copies];
                for (int i = 0; i < copies; ++i) {
                    vector[i] = (Double) gp.nextRandom();
                }
            } else {
                vector = (double[]) gp.nextRandom();
            }

            return new DrawResult(vector, offset);
        }

        private final GaussianProcessRandomGenerator gp;
        private final int copies;
        private final int offset;
        private final boolean isUnivariate;
    }

    private List<Callable<DrawResult>> createTasks() {

        List<Callable<DrawResult>> callers = new ArrayList<Callable<DrawResult>>();

        int offset = 0;
        int index = 0;
        for (GaussianProcessRandomGenerator gp : gpList) {
            final int copies = copyList.get(index);
            if (likelihoodList.get(index) instanceof DistributionLikelihood) { // Univariate
                callers.add(new DrawCaller(gp, copies, offset, true));
                offset += copies;
            } else {
                for (int i = 0; i < copies; ++i) {
                    callers.add(new DrawCaller(gp, 1, offset, false));
                    offset += gp.getDimension();
                }
            }
        }

        return callers;
    }

    @Override
    public Object nextRandom() {
        if (USE_POOL) {
            return nextRandomParallel();
        } else {
            return nextRandomSerial();
        }
    }

    private Object nextRandomParallel() {

        double[] vector = new double[getDimension()];

        try {
            List<Future<DrawResult>> results = pool.invokeAll(callers);

            for (Future<DrawResult> result : results) {
                DrawResult dr = result.get();
                System.arraycopy(dr.result, 0, vector, dr.offset, dr.result.length);
            }

        } catch (InterruptedException e) {
            e.printStackTrace();
        } catch (ExecutionException e) {
            e.printStackTrace();
        }

        return vector;
    }

    private Object nextRandomSerial() {

        int size = 0;
        List<double[]> randomList = new ArrayList<double[]>();
        int index = 0;
        for (GaussianProcessRandomGenerator gp : gpList) {
            final int copies = copyList.get(index);
            if (likelihoodList.get(index) instanceof DistributionLikelihood) { // Univariate
                double[] vector = new double[copies];
                for (int i = 0; i < copies; ++i) {
                    vector[i] = (Double) gp.nextRandom();
                }
                randomList.add(vector);
                size += vector.length;
            } else {
                for (int i = 0; i < copyList.get(index); ++i) {
                    double[] vector = (double[]) gp.nextRandom();
                    randomList.add(vector);
                    size += vector.length;
                }
            }
            ++index;
        }

        double[] result = new double[size];
        int offset = 0;
        for (double[] vector : randomList) {
            System.arraycopy(vector, 0, result, offset, vector.length);
            offset += vector.length;
        }
        return result;
    }

    @Override
    public double logPdf(Object x) {
        throw new RuntimeException("Not yet implemented");
    }
}
