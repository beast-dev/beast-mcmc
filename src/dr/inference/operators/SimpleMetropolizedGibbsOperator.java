/*
 * SimpleMetropolizedGibbsOperator.java
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

/**
 *
 */
package dr.inference.operators;

import dr.inference.model.Likelihood;
import dr.inference.model.Model;

import java.util.logging.Logger;

/**
 * @author Sebastian Hoehna
 */
public abstract class SimpleMetropolizedGibbsOperator extends SimpleOperator implements GeneralOperator {

    /**
     *
     */
    public SimpleMetropolizedGibbsOperator() {
        // Do nothing
    }

    public abstract double doOperation(Likelihood likelihood)
            ;

    /*
      * (non-Javadoc)
      *
      * @see dr.inference.operators.GibbsOperator#getStepCount()
      */
    public abstract int getStepCount();

    /*
      * (non-Javadoc)
      *
      * @see dr.inference.operators.SimpleOperator#getOperatorName()
      */
    @Override
    public abstract String getOperatorName();

    public final double operate() {
        return operate(null);
    }

    public final double operate(Likelihood likelihood) {
        if (operateAllowed) {
            operateAllowed = false;
            return doOperation(likelihood);
        } else
            throw new RuntimeException(
                    "Operate called twice without accept/reject in between!");
    }

    /*
      * (non-Javadoc)
      *
      * @see dr.inference.operators.MCMCOperator#getMaximumAcceptanceLevel()
      */
    public final double getMaximumAcceptanceLevel() {
        return 1.0;
    }

    /*
      * (non-Javadoc)
      *
      * @see dr.inference.operators.MCMCOperator#getMaximumGoodAcceptanceLevel()
      */
    public final double getMaximumGoodAcceptanceLevel() {
        return 1.0;
    }

    /*
      * (non-Javadoc)
      *
      * @see dr.inference.operators.MCMCOperator#getMinimumAcceptanceLevel()
      */
    public final double getMinimumAcceptanceLevel() {
        return 0.005;
    }

    /*
      * (non-Javadoc)
      *
      * @see dr.inference.operators.MCMCOperator#getMinimumGoodAcceptanceLevel()
      */
    public final double getMinimumGoodAcceptanceLevel() {
        return 0.01;
    }

    /*
      * (non-Javadoc)
      *
      * @see dr.inference.operators.MCMCOperator#getPerformanceSuggestion()
      */
    public final String getPerformanceSuggestion() {
        return "";
    }

    /*
      * (non-Javadoc)
      *
      * @see dr.inference.operators.MCMCOperator#getTargetAcceptanceProbability()
      */
    public final double getTargetAcceptanceProbability() {
        return 1.0;
    }

    protected double evaluate(Likelihood likelihood, double pathParameter) {

        double logPosterior = 0.0;

        final double logLikelihood = likelihood.getLogLikelihood() * pathParameter;

        if (Double.isNaN(logLikelihood)) {
            return Double.NEGATIVE_INFINITY;
        }
        // System.err.println("** " + logPosterior + " + " + logLikelihood + " =
        // " + (logPosterior + logLikelihood));
        logPosterior += logLikelihood;

        return logPosterior;
    }

    protected void restore(Likelihood likelihood,
                           Model currentModel, MCMCOperator mcmcOperator, double oldScore) {
        currentModel.restoreModelState();

        // This is a test that the state is correctly restored. The restored
        // state is fully evaluated and the likelihood compared with that before
        // the operation was made.
        likelihood.makeDirty();
        final double testScore = evaluate(likelihood, 1.0);

        if (Math.abs(testScore - oldScore) > 1e-6) {
            Logger.getLogger("error").severe(
                    "State was not correctly restored after reject step.\n"
                            + "Likelihood before: " + oldScore
                            + " Likelihood after: " + testScore + "\n"
                            + "Operator: " + mcmcOperator + " "
                            + mcmcOperator.getOperatorName());
        }
    }

}
