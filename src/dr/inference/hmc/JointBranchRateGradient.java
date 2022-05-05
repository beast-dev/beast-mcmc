/*
 * JointGradient.java
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

import dr.evomodel.treedatalikelihood.continuous.BranchRateGradient;
import dr.evomodel.treedatalikelihood.discrete.BranchRateGradientForDiscreteTrait;
import dr.inference.model.DerivativeOrder;
import dr.inference.model.Likelihood;
import dr.xml.*;

import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.*;

/**
 * @author Alexander Fisher
 */

public class JointBranchRateGradient extends JointGradient {

    private static final boolean COMPUTE_IN_PARALLEL = true;
    private final ExecutorService pool;
    private final List<Callable<double[]>> derivativeCaller;

    private final static String JOINT_BRANCH_RATE_GRADIENT = "JointBranchRateGradient";

    public JointBranchRateGradient(List<GradientWrtParameterProvider> derivativeList) {
        super(derivativeList);

        if (COMPUTE_IN_PARALLEL && derivativeList.size() > 1) {
            pool = Executors.newFixedThreadPool(derivativeList.size());
            derivativeCaller = new ArrayList<>(derivativeList.size());

            for (int i = 0; i < derivativeList.size(); ++i) {
                derivativeCaller.add(new DerivativeCaller(derivativeList.get(i), i));
            }

        } else {
            pool = null;
            derivativeCaller = null;
        }
    }

    @Override
    double[] getDerivativeLogDensity(DerivativeType derivativeType) {

        if (COMPUTE_IN_PARALLEL && pool != null) {
            return getDerivativeLogDensityInParallel(derivativeType);
        } else {
            return super.getDerivativeLogDensity(derivativeType);
        }
    }


    private double[] getDerivativeLogDensityInParallel(DerivativeType derivativeType) {

        if (derivativeType != DerivativeType.GRADIENT) {
            throw new RuntimeException("Not yet implemented");
        }

        final int length = derivativeList.get(0).getDimension();

        double[] derivative = new double[length];

        try {
            List<Future<double[]>> results = pool.invokeAll(derivativeCaller);

            for (Future<double[]> result : results) {
                double[] d = result.get();
                for (int j = 0; j < length; ++j) {
                    derivative[j] += d[j];
                }
            }

        } catch (InterruptedException e) {
            e.printStackTrace();
        } catch (ExecutionException e) {
            e.printStackTrace();
        }

        return derivative;
    }

    // **************************************************************
    // XMLObjectParser
    // **************************************************************

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return JOINT_BRANCH_RATE_GRADIENT;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            List<GradientWrtParameterProvider> derivativeList = new ArrayList<>();

            for (int i = 0; i < xo.getChildCount(); i++) {
                GradientWrtParameterProvider grad = (GradientWrtParameterProvider) xo.getChild(i);
                derivativeList.add(grad);
            }

            return new JointGradient(derivativeList);
        }
        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
                new OrRule(
                        new ElementRule(BranchRateGradient.class, 1, Integer.MAX_VALUE),
                        new ElementRule(BranchRateGradientForDiscreteTrait.class, 1, Integer.MAX_VALUE)
                ),
        };

        public String getParserDescription() {
            return "Joint branch rate gradient";
        }

        public Class getReturnType() {
            return JointBranchRateGradient.class;
        }
    };

    class DerivativeCaller implements Callable<double[]> {

        public DerivativeCaller(GradientWrtParameterProvider gradient, int index) {
            this.gradient = gradient;
            this.index = index;
        }

        public double[] call() throws Exception {
            if (DEBUG_PARALLEL_EVALUATION) {
                System.err.println("Invoking thread #" + index + " for " + gradient.getLikelihood().getId());
            }

            return gradient.getGradientLogDensity();
        }

        private final GradientWrtParameterProvider gradient;
        private final int index;
    }

    public static final boolean DEBUG_PARALLEL_EVALUATION = false;
}
