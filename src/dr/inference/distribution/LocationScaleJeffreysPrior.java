/*
 * LocationScaleJeffreysPrior.java
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

package dr.inference.distribution;

import dr.inference.model.AbstractModelLikelihood;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;
import org.w3c.dom.Document;
import org.w3c.dom.Element;

/**
 * @author Marc A. Suchard
 * @author Robert E. Weiss
 */

public class LocationScaleJeffreysPrior extends AbstractModelLikelihood {


    public LocationScaleJeffreysPrior(Parameter location, Parameter sigma, Parameter gamma,
                                      Parameter alpha0, Parameter beta0, Type type) {

        super("RubioSteele");

        this.location = location;
        this.sigma = sigma;
        this.gamma = gamma;
        this.alpha0 = alpha0;
        this.beta0 = beta0;

        addVariable(location);
        addVariable(sigma);
        addVariable(gamma);

        if (alpha0 != null) addVariable(alpha0);
        if (beta0 != null) addVariable(beta0);

        this.type = type;

    }

    // **************************************************************
    // Likelihood IMPLEMENTATION
    // **************************************************************

    public Model getModel() {
        return this;
    }

    /**
     * Calculate the log likelihood of the current state.
     *
     * @return the log likelihood.
     */
    public double getLogLikelihood() {
        double a0 = alpha0 != null ? alpha0.getParameterValue(0) : 0.0;
        double b0 = beta0 != null ? beta0.getParameterValue(0) : 0.0;
        return type.logPdf(location.getParameterValue(0), sigma.getParameterValue(0), gamma.getParameterValue(0),
                a0, b0);
    }

    public void makeDirty() {
    }

    public void acceptState() {
        // DO NOTHING
    }

    public void restoreState() {
        // DO NOTHING
    }

    public void storeState() {
        // DO NOTHING
    }

    protected void handleModelChangedEvent(Model model, Object object, int index) {
        // DO NOTHING
    }

    protected final void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        // DO NOTHING
    }

    // **************************************************************
    // XMLElement IMPLEMENTATION
    // **************************************************************

    public Element createElement(Document d) {
        throw new RuntimeException("Not implemented yet!");
    }

    private final Parameter location;
    private final Parameter sigma;
    private final Parameter gamma;
    private final Parameter alpha0;
    private final Parameter beta0;

    private final Type type;

    private interface Compute {
        double logPdf(double mu, double sigma, double gamma, double alpha0, double beta0);
    }

    public enum Type implements Compute {

        ISF_AG_BETA("inverseScaleFactorsAGBeta", new Compute() {
            public double logPdf(double mu, double sigma, double gamma, double alpha0, double beta0) {
                return -Math.log(sigma) + (2 * alpha0 - 1) * Math.log(gamma) - (alpha0 + beta0) * Math.log(1.0 + gamma * gamma);
            }
        }),

        EPSILON_SKEW_AG_BETA("epsilonSkewAGBeta", new Compute() {
            public double logPdf(double mu, double sigma, double gamma, double alpha0, double beta0) {
                return -Math.log(sigma) + (beta0 - 1) * Math.log(1.0 - gamma) + (alpha0 - 1) * Math.log(1.0 + gamma);
            }
        }),

        LOGISTIC_AG_BETA("logisticAGBeta", new Compute() {
            public double logPdf(double mu, double sigma, double gamma, double alpha0, double beta0) {
                return -Math.log(sigma)
                        + alpha0 * Math.log(1.0 + Math.exp(2.0 * gamma))
                        + beta0 * Math.log(1.0 + Math.exp(-2.0 * gamma))
                        - (alpha0 + beta0) * Math.log(1.0 + Math.cosh(2.0 * gamma));
            }
        }),

        TWO_SCALE("twoScale", new Compute() {
            public double logPdf(double mu, double sigma, double gamma, double alpha0, double beta0) {
                return -(Math.log(sigma) + Math.log(gamma) + Math.log(sigma + gamma));
            }
        }),

        TWO_SCALE_INDEPENDENT("twoScaleIndependent", new Compute() {
            public double logPdf(double mu, double sigma, double gamma, double alpha0, double beta0) {
                throw new RuntimeException("Not yet implemented");
            }
        }),

        EPSILON_SKEW("epsilonSkew", new Compute() {
            @Override
            public double logPdf(double mu, double sigma, double gamma, double alpha0, double beta0) {
                return -2.0 * Math.log(sigma) - Math.log(1.0 - gamma * gamma);
            }
        }),

        EPSILON_SKEW_INDEPENDENT("epsilonSkewIndependent", new Compute() {
            @Override
            public double logPdf(double mu, double sigma, double gamma, double alpha0, double beta0) {
                return -Math.log(sigma) - 0.5 * Math.log(1.0 - gamma * gamma);
            }
        }),

        TWO_SCALE2("twoScale", new Compute() {
            public double logPdf(double mu, double sigma, double gamma, double alpha0, double beta0) {
                return -(Math.log(sigma) + Math.log(gamma) + Math.log(sigma + gamma));
            }
        });

        Type(String text, Compute compute) {
            this.text = text;
            this.compute = compute;
        }

        private final String text;
        private final Compute compute;

        public static Type parseFromString(String text) {
            for (Type scheme : Type.values()) {
                if (scheme.toString().compareToIgnoreCase(text) == 0)
                    return scheme;
            }
            return null;
        }

        public String toString() {
            return text;
        }

        public double logPdf(double mu, double sigma, double gamma, double alpha0, double beta0) {
            return compute.logPdf(mu, sigma, gamma, alpha0, beta0);
        }
    }
}

