/*
 * TwoPieceLocationScaleDistributionModel.java
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

import dr.inference.model.AbstractModel;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;
import dr.inferencexml.distribution.TwoPieceLocationScaleDistributionModelParser;
import dr.math.UnivariateFunction;
import dr.math.distributions.Distribution;
import dr.math.distributions.RandomGenerator;
import org.w3c.dom.Document;
import org.w3c.dom.Element;

/**
 * @author Marc A. Suchard
 * @author Robert E. Weiss
 */

public class TwoPieceLocationScaleDistributionModel extends AbstractModel implements ParametricDistributionModel, RandomGenerator {

    public TwoPieceLocationScaleDistributionModel(Parameter locationParam, Distribution distribution,
                                                  Parameter sigmaParameter,
                                                  Parameter gammaParameter, Parameterization parameterization) {
        super(TwoPieceLocationScaleDistributionModelParser.DISTRIBUTION_MODEL);
        this.locationParameter = locationParam;
        this.sigmaParameter = sigmaParameter;
        this.gammaParameter = gammaParameter;
        this.distribution = distribution;

        addVariable(locationParam);
        addVariable(sigmaParameter);
        addVariable(gammaParameter);

        locationParam.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, Double.NEGATIVE_INFINITY, 1));
        sigmaParameter.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, 1));
        gammaParameter.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, 1));

        this.parameterization = parameterization;

        // TODO Upgrade Distribution to DistributionModel
    }

    private double getLowerScale() {
        return parameterization.getLowerScale(sigmaParameter.getParameterValue(0), gammaParameter.getParameterValue(0));
    }

    private double getUpperScale() {
        return parameterization.getUpperScale(sigmaParameter.getParameterValue(0), gammaParameter.getParameterValue(0));
    }

    private double getLocation() {
        return locationParameter.getParameterValue(0);
    }

    // *****************************************************************
    // Interface Distribution
    // *****************************************************************

    public double pdf(double x) {
        return Math.exp(logPdf(x));
    }

    public double logPdf(double x) {
        final double mu = getLocation();
        final double lowerScale = getLowerScale();
        final double upperScale = getUpperScale();

        double t = x - mu;
        if (x < mu) {
            t /= lowerScale;
        } else {
            t /= upperScale;
        }

        return Math.log(2.0) - Math.log(lowerScale + upperScale) + distribution.logPdf(t);
    }

    public double cdf(double x) {
        throw new RuntimeException("Not yet implemented.");
    }

    public double quantile(double y) {
        throw new RuntimeException("Not yet implemented.");
    }

    public double mean() {
        throw new RuntimeException("Not yet implemented.");
    }

    public double variance() {
        throw new RuntimeException("Not yet implemented.");
    }

    public final UnivariateFunction getProbabilityDensityFunction() {
        throw new RuntimeException("Not yet implemented.");
    }

    public Object nextRandom() {
        throw new RuntimeException("Not yet implemented.");
    }

    public double logPdf(Object x) {
        double v = (Double) x;
        return logPdf(v);
    }

    // *****************************************************************
    // Interface DensityModel
    // *****************************************************************

    @Override
    public double logPdf(double[] x) {
        return logPdf(x[0]);
    }

    @Override
    public Variable<Double> getLocationVariable() {
        throw new UnsupportedOperationException("Not implemented");
    }

    // *****************************************************************
    // Interface Model
    // *****************************************************************

    public void handleModelChangedEvent(Model model, Object object, int index) {
        // no intermediates need to be recalculated...
    }

    protected final void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        // no intermediates need to be recalculated...
    }

    protected void storeState() {
    } // no additional state needs storing

    protected void restoreState() {
    } // no additional state needs restoring

    protected void acceptState() {
    } // no additional state needs accepting

    public Element createElement(Document document) {
        throw new RuntimeException("Not yet implemented!");
    }

    // **************************************************************
    // Private instance variables
    // **************************************************************

    private final Parameter locationParameter;
    private final Parameter sigmaParameter;
    private final Parameter gammaParameter;

    private final Distribution distribution;
    private final Parameterization parameterization;


    private interface Scale {
        double getLowerScale(double sigma, double gamma);

        double getUpperScale(double sigma, double gamma);
    }

    public enum Parameterization implements Scale {

        TWO_SCALE("twoScale", new Scale() {
            public double getLowerScale(double sigma, double gamma) {
                return sigma;
            }

            public double getUpperScale(double sigma, double gamma) {
                return gamma;
            }
        }),

        ISF("inverseScaleFactors", new Scale() {
            public double getLowerScale(double sigma, double gamma) {
                return sigma / gamma;
            }

            public double getUpperScale(double sigma, double gamma) {
                return sigma * gamma;
            }
        }),

        EPSILON_SKEW("epsilonSkew", new Scale() {
            public double getLowerScale(double sigma, double gamma) {
                return sigma * (1.0 - gamma);
            }

            public double getUpperScale(double sigma, double gamma) {
                return sigma * (1.0 + gamma);
            }
        }),

        LOGISTIC_AG("logisticAG", new Scale() {
            public double getLowerScale(double sigma, double gamma) {
                return sigma * (1.0 + Math.exp(-2.0 * gamma));
            }

            public double getUpperScale(double sigma, double gamma) {
                return sigma * (1.0 + Math.exp(2.0 * gamma));
            }
        });

        Parameterization(String text, Scale scale) {
            this.text = text;
            this.scale = scale;
        }

        private final String text;
        private final Scale scale;

        public static Parameterization parseFromString(String text) {
            for (Parameterization scheme : Parameterization.values()) {
                if (scheme.toString().compareToIgnoreCase(text) == 0)
                    return scheme;
            }
            return null;
        }

        public String toString() {
            return text;
        }

        public double getLowerScale(double sigma, double gamma) {
            return scale.getLowerScale(sigma, gamma);
        }

        public double getUpperScale(double sigma, double gamma) {
            return scale.getUpperScale(sigma, gamma);
        }

    }

}
