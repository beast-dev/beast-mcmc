/*
 * GaussianProcessKernel.java
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

package dr.math.distributions.gp;

import dr.inference.model.*;

import java.util.List;

/**
 * @author Marc A. Suchard
 * @author Filippo Monti
 */
public interface GaussianProcessKernel {

    double getUnscaledCovariance(double x, double y);

    @SuppressWarnings("unused")
    double getUnscaledCovariance(double[] x, double[] y);

    double getScale();

    List<Parameter> getParameters();

    double getUnscaledFirstDerivative(double x, double y);
    double getUnscaledSecondDerivative(double x, double y);
    double getNormalizationFactor(double[] x);
    default boolean isUnitaryVariance() { return false; }

    class KernelDerivatives extends Base {
        GaussianProcessKernel kernel;
        boolean doSecondDerivative = false;
        public KernelDerivatives(GaussianProcessKernel kernel, boolean doSecondDerivative) {
            super("KernelFirstDerivative", kernel.getParameters(), kernel.isUnitaryVariance());
            this.kernel = kernel;
            this.doSecondDerivative = doSecondDerivative;
        }
        @Override
        public double getUnscaledCovariance(double x, double y) {
            if (!doSecondDerivative) {
                return kernel.getUnscaledFirstDerivative(x, y);
            } else {
                return kernel.getUnscaledSecondDerivative(x, y);
            }
        }

        @Override
        public double getUnscaledCovariance(double[] x, double[] y) {
            throw new RuntimeException("Method not yet implemented");
        }
    }

    class Linear extends Base {

        public Linear(String name, List<Parameter> parameters, boolean unitaryVariance) {
            super(name, parameters, unitaryVariance);
        }

        public double getUnscaledCovariance(double x, double y) {
            return  x * y;
        }

        public double getUnscaledCovariance(double[] x, double[] y) {
            final int dim = x.length;

            double product = 0.0;
            for (int i = 0; i < dim; ++i) {
                product += x[i] * y[i];
            }

            return product;
        }

        private static final String TYPE = "DotProduct";

        public double getUnscaledFirstDerivative(double x, double y) {
            return y;
        }
        public double getUnscaledSecondDerivative(double x, double y) {
           return 1;
        }

        @Override
        public double getNormalizationFactor(double[] x) {
            double normalizationConstant = 0.0;
            for (int i = 0; i < x.length; i++) {
                normalizationConstant += x[i] * x[i];
            }
            return normalizationConstant / x.length;
        }

        @Override
        public double computeGradientWrtLength(double a, double b, double l) {
            throw new RuntimeException("The linear kernel does not have a length parameter");
        }

    }
    class NeuralNet extends Base {
        private static final String TYPE = "NeuralNet";
        public NeuralNet(String name, List<Parameter> parameters, boolean unitaryVariance) {
            super(name, parameters, unitaryVariance);
        }

        @Override
        public double getUnscaledCovariance(double x, double y) {
            return getUnscaledCovariance(new double[]{x}, new double[]{y});
        }

        @Override
        public double getUnscaledCovariance(double[] x, double[] y) {
            if (x.length != y.length) {
                throw new IllegalArgumentException("Input vectors must have the same length.");
            }

            double dot = 0.0;
            double normX = 0.0;
            double normY = 0.0;

            for (int i = 0; i < x.length; i++) {
                dot += x[i] * y[i];
                normX += x[i] * x[i];
                normY += y[i] * y[i];
            }

            double numerator = 2 * dot;
            double denominator = Math.sqrt((1 + 2 * normX) * (1 + 2 * normY));
            double argument = numerator / denominator;

            // Clamp to [-1, 1] for numerical safety
            argument = Math.max(-1.0, Math.min(1.0, argument));

            return (2.0 / Math.PI) * Math.asin(argument);
        }

        @Override
        public double computeGradientWrtLength(double a, double b, double l) {
            throw new RuntimeException("The linear kernel does not have a length parameter");
        }
    }

    class SquaredExponential extends Base {

        public SquaredExponential(String name, List<Parameter> parameters, boolean unitaryVariance) {
            super(name, parameters, unitaryVariance);
        }

        private double functionalForm(double normSquared) {
            double length = getLength();
            return Math.exp(-normSquared / (2 * length * length));
        }

        public double getUnscaledCovariance(double x, double y) {
            double diff = x - y;
            return functionalForm(diff * diff);
        }

        public double getUnscaledCovariance(double[] x, double[] y) {
            return functionalForm(normSquared(x, y));
        }

        private static final String TYPE = "SquaredExponential";

        public double getUnscaledFirstDerivative(double x, double y) { //wrt x
            double diff = x - y;
            double length = getLength();
            return - diff / (length * length) * getUnscaledCovariance(x, y);
        }
        public double getUnscaledSecondDerivative(double x, double y) {
            double diff = x - y;
            double norm = diff * diff;
            double length = getLength();
            double lengthSquared = length * length;
            return (1/lengthSquared -  norm / (lengthSquared * lengthSquared)) * getUnscaledCovariance(x, y);
        }

        @Override
        public double computeGradientWrtLength(double a, double b, double l) {
            double normSquared = (a - b) * (a - b);
            double scale = getScale();
            return scale * Math.exp(-normSquared / (2 * l * l)) *
                    normSquared / (l * l * l);
        }
    }

    class OrnsteinUhlenbeck extends NormedBase {

        public OrnsteinUhlenbeck(String name, List<Parameter> parameters, boolean unitaryVariance) { 
            super(name, parameters, unitaryVariance); 
        }

        double functionalForm(double norm) {
            double length = getLength();
            return Math.exp(-norm / length);
        }
        
        private static final String TYPE = "OrnsteinUhlenbeck";
    }

    class MaternFiveHalves extends NormedBase {

        public MaternFiveHalves(String name, List<Parameter> parameters, boolean unitaryVariance) { 
            super(name, parameters, unitaryVariance); 
        }

        double functionalForm(double norm) {
            double length = getLength();

            double argument1 = Math.sqrt(5) * norm / length;
            double argument2 = 5 * norm * norm / (3 * length * length);

            return (1 + argument1 + argument2) * Math.exp(-argument1);
        }

        private static final String TYPE = "Matern5/2";
    }

    class MaternThreeHalves extends NormedBase {

        public MaternThreeHalves(String name, List<Parameter> parameters, boolean unitaryVariance) { 
            super(name, parameters, unitaryVariance); }

        double functionalForm(double norm) {
            double length = getLength();

            double argument = Math.sqrt(3) * norm / length;

            return (1 + argument) * Math.exp(-argument);
        }

        private static final String TYPE = "Matern3/2";
    }

    abstract class NormedBase extends Base {

        public NormedBase(String name, List<Parameter> parameters, boolean unitaryVariance) { 
            super(name, parameters, unitaryVariance); 
        }

        abstract double functionalForm(double norm);

        public double getUnscaledCovariance(double x, double y) {
            double norm = Math.abs(x - y);
            return functionalForm(norm);
        }

        public double getUnscaledCovariance(double[] x, double[] y) {
            double norm = Math.sqrt(normSquared(x, y));
            return functionalForm(norm);
        }
    }
    
    enum AllKernels {
        LINEAR(Linear.TYPE) {
            GaussianProcessKernel factory(String name, List<Parameter> parameters, boolean unitaryVariance) {
                return new Linear(name, parameters, unitaryVariance); 
            }
        },
        NEURAL_NET(NeuralNet.TYPE) {
            GaussianProcessKernel factory(String name, List<Parameter> parameters, boolean unitaryVariance) {
                return new NeuralNet(name, parameters, unitaryVariance); 
            }
        },
        OU(OrnsteinUhlenbeck.TYPE) {
            GaussianProcessKernel factory(String name, List<Parameter> parameters, boolean unitaryVariance) {
                return new OrnsteinUhlenbeck(name, parameters, unitaryVariance); 
            }
        },
        SE(SquaredExponential.TYPE) {
            GaussianProcessKernel factory(String name, List<Parameter> parameters, boolean unitaryVariance) {
                return new SquaredExponential(name, parameters, unitaryVariance); 
            }
        },
        MATERN32(MaternThreeHalves.TYPE) {
            GaussianProcessKernel factory(String name, List<Parameter> parameters, boolean unitaryVariance) {
                return new MaternThreeHalves(name, parameters, unitaryVariance); 
            }
        },
        MATERN52(MaternFiveHalves.TYPE) {
            GaussianProcessKernel factory(String name, List<Parameter> parameters, boolean unitaryVariance) {
                return new MaternFiveHalves(name, parameters, unitaryVariance); 
            }
        };

        AllKernels(String name) { this.name = name; }

        abstract GaussianProcessKernel factory(String name, List<Parameter> parameters, boolean unitaryVariance); 
        public String toString() { return name; }
        private final String name;
    }
    static GaussianProcessKernel factory(String type, String name, List<Parameter> parameters, boolean unitaryVariance)
            throws IllegalArgumentException {
        for (AllKernels kernel : AllKernels.values()) {
            if (type.equalsIgnoreCase(kernel.toString())) {
                return kernel.factory(name, parameters, unitaryVariance); 
            }
        }
        throw new IllegalArgumentException("Unknown kernel type");
    }

    @FunctionalInterface
    interface HyperparameterGradientFunction {
        double apply(double a, double b, double hyperValue);
    }

    abstract class Base extends AbstractModel implements GaussianProcessKernel {

        final List<Parameter> parameters;
        final boolean unitaryVariance;

        public Base(String name,
                    List<Parameter> parameters,
                    boolean unitaryVariance) {
            super(name);
            this.parameters = parameters;
            this.unitaryVariance = unitaryVariance;

            for (Parameter p : parameters) {
                addVariable(p);
            }
        }

        double normSquared(double[] x, double[] y) {
            final int dim = x.length;

            double normSquared = 0.0;
            for (int i = 0; i < dim; ++i) {
                double diff = x[i] - y[i];
                normSquared += diff * diff;
            }

            return normSquared;
        }
        // this is a normalization factor necessary to get variance 1  for non-stationary kernels
        public double getNormalizationFactor(double[] x) {
            return 1.0;
        }
        
        @Override
        public boolean isUnitaryVariance() {
            return unitaryVariance;
        }

        @Override
        public double getScale() {
            return parameters.get(0).getParameterValue(0);
        }

        double getLength() {return parameters.get(1).getParameterValue(0);}

        public List<Parameter> getParameters() {
            return parameters;
        }

        public HyperparameterGradientFunction getScaleGradientFunction() {
            return (a, b, hyperValue) -> getUnscaledCovariance(a, b);
        }
        public HyperparameterGradientFunction getLengthGradientFunction() {
            return this::computeGradientWrtLength;
        }

        public double getUnscaledFirstDerivative(double x, double y) {
            throw new RuntimeException("Method implemented in the subclasses");
        }
        public double getUnscaledSecondDerivative(double x, double y) {
            throw new RuntimeException("Method implemented in the subclasses");
        }

        public double computeGradientWrtLength(double a, double b, double hyperValue) {
            throw new RuntimeException("Method implemented in the subclasses");
        }

        @Override
        protected void handleModelChangedEvent(Model model, Object object, int index) { }

        @Override
        protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
            if (parameters.contains((Parameter) variable)) {
                fireModelChanged(variable, index);
            } else {
                throw new IllegalArgumentException("Unknown variable");
            }
        }

        @Override
        protected void storeState() { }

        @Override
        protected void restoreState() { }

        @Override
        protected void acceptState() { }
    }
}
