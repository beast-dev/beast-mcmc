/*
 * NormalGammaDistribution.java
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


import dr.math.functionEval.GammaFunction;

/**
 * Normal-gamma distribution
 * @author Matthew Hall
 */
public class NormalGammaDistribution implements MultivariateDistribution {

    public final String TYPE = "normalGammaDistribution";

    private double mu;
    private double lambda;
    private double alpha;
    private double beta;

    public NormalGammaDistribution(double mu, double lambda, double alpha, double beta){
        this.mu=mu;
        this.lambda=lambda;
        this.alpha=alpha;
        this.beta=beta;
    }

    public double pdf(double[] x){
        return Math.pow(beta, alpha)*Math.sqrt(lambda)*Math.pow(x[1],alpha-0.5)*Math.exp(-beta*x[1])
                *Math.exp(-lambda*x[1]*Math.pow(x[0]-mu,2)/2)/(GammaFunction.gamma(alpha)*Math.sqrt(2*Math.PI));
    }

    public double logPdf(double[] x) {
        return alpha*Math.log(beta) + 0.5*Math.log(lambda) + (alpha-0.5)*Math.log(x[1])-beta*x[1]
                -lambda*x[1]*Math.pow(x[0]-mu,2)/2 - GammaFunction.logGamma(alpha) - 0.5*Math.log(2*Math.PI);
    }

    public double[][] getScaleMatrix() {
        throw new RuntimeException("Not implemented");
    }

    public double[] getMean() {
        return new double[]{mu, alpha/beta};
    }

    public String getType() {
        return TYPE;
    }

    public double[] getParameters() {
        return new double[]{mu, lambda, alpha, beta};
    }
}
