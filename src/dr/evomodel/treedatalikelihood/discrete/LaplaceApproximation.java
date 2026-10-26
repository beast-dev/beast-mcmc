/*
 * LaplaceApproximation.java
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

package dr.evomodel.treedatalikelihood.discrete;

import com.github.lbfgs4j.liblbfgs.Function;
import dr.math.distributions.MultivariateNormalDistribution;
import dr.util.Timer;
import dr.xml.Reportable;

import java.util.List;

/**
 * A Gaussian approximation of the maximizer's objective at its optimum, in the coordinates it optimizes (y): the mode,
 * the covariance A^-1 and, when the maximizer includes the Jacobian, the Laplace log evidence. The curvature and the
 * evidence come from LaplaceEvidence (plan 8 X3), which corrects the earlier Hessian taken in the untransformed
 * parameter and the sign of the log determinant.
 *
 * @author Andy Magee
 */

public class LaplaceApproximation implements Reportable {

    private final MaximizerWrtParameter maximizer;
    private final boolean diagonal;
    private final boolean estimateKL;
    private final boolean computeML;
    private LaplaceEvidence laplace;
    private double time = 0.0;
    private double KL = Double.NaN;
    private double ML = Double.NaN;

    public LaplaceApproximation(MaximizerWrtParameter maximizer, boolean diagonal, boolean KL, boolean ML) {
        this.maximizer = maximizer;
        this.diagonal = diagonal;
        this.estimateKL = KL;
        this.computeML = ML;
    }

    public void approximate() {
        Timer timer = new Timer();
        timer.start();
        laplace = new LaplaceEvidence(maximizer);
        timer.stop();
        time += timer.toSeconds();

        if ( estimateKL ) {
            // KL(approximation || objective), up to the objective's normalizing constant, over the free coordinates
            int[] coordinates = laplace.getCoordinates();
            double[] mode = laplace.getMode();
            double[] mean = new double[coordinates.length];
            for (int i = 0; i < coordinates.length; i++) {
                mean[i] = mode[coordinates[i]];
            }
            MultivariateNormalDistribution distribution = new MultivariateNormalDistribution(mean, laplace.getPrecision());
            Function function = maximizer.getFunction();
            int nSamples = 1000;
            KL = 0.0;
            for (int i = 0; i < nSamples; i++) {
                double[] draw = distribution.nextMultivariateNormal();
                double[] y = mode.clone();
                for (int j = 0; j < coordinates.length; j++) {
                    y[coordinates[j]] = draw[j];
                }
                double lnPosterior = -function.valueAt(y);
                KL += distribution.logPdf(draw) - lnPosterior;
            }
            KL /= nSamples;
            function.valueAt(mode);
        }

        if ( computeML ) {
            ML = laplace.getLogEvidence();
        }
    }

    /**
     * @return the Laplace log evidence (computeML), or NaN
     */
    public double getLogMarginalLikelihood() {
        return ML;
    }

    public LaplaceEvidence getLaplaceEvidence() {
        return laplace;
    }

    @Override
    public String getReport() {
        StringBuilder sb = new StringBuilder();

        if ( !maximizer.wasExecuted() || laplace == null ) {
            sb.append("Not yet executed.");
        } else {
            sb.append("Laplace approximation in the maximizer's transformed coordinates (y): the mode and the ")
                    .append(diagonal ? "standard deviations" : "covariance").append(" are of y, not of the parameter.\n");
            int[] coordinates = laplace.getCoordinates();
            double[] mode = laplace.getMode();
            double[][] covariance = laplace.getCovariance();
            List<MaximizerWrtParameter.Coordinate> optimum = maximizer.getOptimum();
            if (coordinates.length < mode.length) {
                sb.append("Coordinates held at a bound or dependent (the simplex's last) are left out, and the ")
                        .append("approximation is conditional on them.\n");
            }
            for (int i = 0; i < coordinates.length; i++) {
                sb.append(optimum.get(coordinates[i]).name).append(": y = ").append(mode[coordinates[i]]);
                if ( diagonal ) {
                    sb.append(", sd(y) = ").append(Math.sqrt(covariance[i][i]));
                }
                sb.append("\n");
            }
            if ( !diagonal ) {
                sb.append("Covariance matrix of y: ").append(new dr.math.matrixAlgebra.Matrix(covariance)).append("\n");
            }
            if (estimateKL) {
                sb.append("KL(approximation||objective), up to its normalizing constant = ").append(KL).append("\n");
            }
            if (computeML) {
                sb.append("log(marginal likelihood), Laplace = ").append(ML).append("\n");
            }
            sb.append("Time to compute Hessian: ").append(time).append("s\n");
        }

        return sb.toString();

    }
}
