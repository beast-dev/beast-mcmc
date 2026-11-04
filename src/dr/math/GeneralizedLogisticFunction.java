/*
 * GeneralizedLogisticFunction.java
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

package dr.math;

import org.apache.commons.math.ConvergenceException;
import org.apache.commons.math.FunctionEvaluationException;
import org.apache.commons.math.analysis.UnivariateRealFunction;
import org.apache.commons.math.analysis.integration.RombergIntegrator;
import org.apache.commons.math.analysis.integration.UnivariateRealIntegrator;

/**
 * @author Marc A. Suchard
 */
public class GeneralizedLogisticFunction {

    public static double evaluate(double x,
                                  double leftTime,
                                  double rightTime,
                                  double leftAsymptote,
                                  double rightAsymptote,
                                  double growthRate,
                                  double inflectionTime,
                                  double inflectionHeight) {
        return evaluate(x, leftTime, rightTime, leftAsymptote, rightAsymptote,
                growthRate, inflectionTime, inflectionHeight, 1.0, 1.0);
    }

    /**
     *  Based on https://en.wikipedia.org/wiki/Generalised_logistic_function
     */
    public static double evaluate(double x,
                                  double leftTime,
                                  double rightTime,
                                  double leftAsymptote, // A
                                  double rightAsymptote, // K
                                  double growthRate, // B
                                  double inflectionTime, // M
                                  double inflectionHeight, // Q
                                  double nu,
                                  double C) {

        final double value;
        if (x <= leftTime) {
            value =  leftAsymptote;
        } else if (x >= rightTime) {
            value = rightAsymptote;
        } else {
            double timeOdds = (rightTime - x) / (x - leftTime);
            double inflectionOdds = (rightTime - inflectionTime) / (inflectionTime - leftTime);
            double inflation = Math.exp(-growthRate * (timeOdds - inflectionOdds));
            double weight = C + inflectionHeight * Math.pow(timeOdds / inflectionOdds, growthRate); //inflation;
            if (nu != 1.0) {
                weight = Math.pow(weight, 1.0 / nu);
            }
            value = leftAsymptote + (rightAsymptote - leftAsymptote) / weight;
        }
        return value;
    }

    public static double integrate(double x0, double x1,
                                   double leftTime,
                                   double rightTime,
                                   double leftAsymptote,
                                   double rightAsymptote,
                                   double growthRate,
                                   double inflectionTime,
                                   double inflectionHeight) {
        return integrate(x0, x1, leftTime, rightTime, leftAsymptote, rightAsymptote,
                growthRate, inflectionTime, inflectionHeight, 1.0, 1.0);
    }

    public static double integrate(double x0, double x1,
                                   double leftTime,
                                   double rightTime,
                                   double leftAsymptote,
                                   double rightAsymptote,
                                   double growthRate,
                                   double inflectionTime,
                                   double inflectionHeight,
                                   double nu,
                                   double C) {

        UnivariateRealFunction f = v -> evaluate(v, leftTime, rightTime,
                leftAsymptote, rightAsymptote, growthRate, inflectionTime, inflectionHeight, nu, C);

        double result = 0.0;
        try {
            result = integrator.integrate(f, x0, x1);
        } catch (ConvergenceException | FunctionEvaluationException e) {
            e.printStackTrace();
        }
        return result;
    }

    private final static UnivariateRealIntegrator integrator = new RombergIntegrator();
}
