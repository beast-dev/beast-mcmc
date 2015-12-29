/*
 * GeneralisedGaussLaguerreQuadrature.java
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

package dr.math;

/**
 * Approximate the integral from min to +INF of x^alpha * e^-Bx * f(x) by generalised Gauss-Riemann quadrature
 *
 * Adapted from Numerical Recipes
 */
public class GeneralisedGaussLaguerreQuadrature implements Integral {

    private double B;
    private double alpha;
    final private int noPoints;
    final private double[] abscissae;
    final private double[] coefficients;

    public GeneralisedGaussLaguerreQuadrature(double B, double alpha, int noPoints){
        this.B = B;
        this.alpha = alpha;
        this.noPoints = noPoints;
        abscissae = new double[noPoints];
        coefficients = new double[noPoints];
        setupArrays();
    }

    public GeneralisedGaussLaguerreQuadrature(int noPoints){
        this.noPoints = noPoints;
        abscissae = new double[noPoints];
        coefficients = new double[noPoints];
    }

    public void setB(double B){
        this.B = B;
        setupArrays();
    }

    public void setAlpha(double alpha){
        this.alpha = alpha;
        setupArrays();
    }

    public void setAlphaAndB(double alpha, double B){
        this.alpha = alpha;
        this.B = B;
        setupArrays();
    }

    private void setupArrays(){
        final int maxIterations = 10;
        final double eps = 1E-14;

        double z = 0;

        double p1=0, p2=0, p3=0, pp=0;

        for(int i=0; i<noPoints; i++){
            if(i==0){
                z = (1.0+alpha)*(3.0+0.92*alpha)/(1.0+2.4*noPoints+1.8*alpha);
            } else if (i==1){
                z += (15.0+6.25*alpha)/(1.0+0.9*alpha+2.5*noPoints);
            } else {
                double ai = i-1;
                z += ((1.0+2.55*ai)/(1.9*ai) + 1.26*ai*alpha/(1.0+3.5*ai))*(z-abscissae[i-2])/(1.0+0.3*alpha);
            }

            boolean failed = true;

            for(int its = 0 ; its<maxIterations; its++){
                p1 = 1.0;
                p2 = 0.0;
                for(int j=0; j<noPoints; j++){
                    p3 = p2;
                    p2 = p1;
                    p1 = ((2*j+1+alpha-z)*p2-(j+alpha)*p3)/(j+1);
                }
                pp = (noPoints*p1-(noPoints+alpha)*p2)/z;
                double z1=z;
                z=z1-p1/pp;
                if(Math.abs(z-z1) <= eps){
                    failed = false;
                    break;
                }
            }
            if(failed){
                throw new RuntimeException("Too many iterations");
            }
            abscissae[i] = z;

            coefficients[i] = -Math.exp(GammaFunction.lnGamma(alpha+noPoints) - GammaFunction.lnGamma((double)noPoints))/
                    (pp*noPoints*p2);

        }

    }


    public double integrate(UnivariateFunction f, double min, double max) {
        if(max!=Double.POSITIVE_INFINITY){
            throw new RuntimeException("Gauss-Laguerre quadrature is for integration with an infinite upper limit");
        }
        else{
            return integrate(f, min);
        }
    }

    public double logIntegrate(UnivariateFunction f, double min, double max) {
        if(max!=Double.POSITIVE_INFINITY){
            throw new RuntimeException("Gauss-Laguerre quadrature is for integration with an infinite upper limit");
        }
        else{
            return logIntegrate(f, min);
        }
    }

    public double integrate(UnivariateFunction f, double min){
        double integral = 0;

        for(int i=0; i<noPoints; i++){
            integral += coefficients[i]*f.evaluate(abscissae[i]/B + min);
        }

        integral *= 1/Math.pow(B, alpha+1);

        return integral;
    }

    public double logIntegrate(UnivariateFunction f, double min){
        try {

            double logIntegral = Double.NEGATIVE_INFINITY;

            for (int i = 0; i < noPoints; i++) {
                logIntegral = LogTricks.logSum(logIntegral, Math.log(coefficients[i]) +
                        ((UnivariateFunction.AbstractLogEvaluatableUnivariateFunction)f)
                                .logEvaluate(min + abscissae[i] / B));
            }

            logIntegral += -(alpha + 1) * Math.log(B);

            return logIntegral;
        } catch(ClassCastException e){
            throw new RuntimeException(f.toString()+" has no logEvaluate method");
        }

    }


}
