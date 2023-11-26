package dr.math.distributions.gp;

import dr.inference.model.AbstractModel;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;

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

    class Linear extends Base {

        public Linear(String name, List<Parameter> parameters) {
            super(name, parameters);
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

    }

    class SquaredExponential extends Base {

        public SquaredExponential(String name, List<Parameter> parameters) {
            super(name, parameters);
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
    }

    class OrnsteinUhlenbeck extends NormedBase {

        public OrnsteinUhlenbeck(String name, List<Parameter> parameters) { super(name, parameters); }

        double functionalForm(double norm) {
            double length = getLength();
            return Math.exp(-norm / length);
        }

        private static final String TYPE = "OrnsteinUhlenbeck";
    }

    class MaternFiveHalves extends NormedBase {

        public MaternFiveHalves(String name, List<Parameter> parameters) { super(name, parameters); }

        double functionalForm(double norm) {
            double length = getLength();

            double argument1 = Math.sqrt(5) * norm / length;
            double argument2 = 5 * norm * norm / (3 * length * length);

            return (1 + argument1 + argument2) * Math.exp(-argument1);
        }

        private static final String TYPE = "Matern5/2";
    }

    class MaternThreeHalves extends NormedBase {

        public MaternThreeHalves(String name, List<Parameter> parameters) { super(name, parameters); }

        double functionalForm(double norm) {
            double length = getLength();

            double argument = Math.sqrt(3) * norm / length;

            return (1 + argument) * Math.exp(-argument);
        }

        private static final String TYPE = "Matern3/2";
    }

    abstract class NormedBase extends Base {

        public NormedBase(String name, List<Parameter> parameters) { super(name, parameters); }

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
            GaussianProcessKernel factory(String name, List<Parameter> parameters) {
                return new Linear(name, parameters);
            }
        },
        OU(OrnsteinUhlenbeck.TYPE) {
            GaussianProcessKernel factory(String name, List<Parameter> parameters) {
                return new OrnsteinUhlenbeck(name, parameters);
            }
        },
        SE(SquaredExponential.TYPE) {
            GaussianProcessKernel factory(String name, List<Parameter> parameters) {
                return new SquaredExponential(name, parameters);
            }
        },
        MATERN32(MaternThreeHalves.TYPE) {
            GaussianProcessKernel factory(String name, List<Parameter> parameters) {
                return new MaternThreeHalves(name, parameters);
            }
        },
        MATERN52(MaternFiveHalves.TYPE) {
            GaussianProcessKernel factory(String name, List<Parameter> parameters) {
                return new MaternFiveHalves(name, parameters);
            }
        };

        AllKernels(String name) { this.name = name; }
        abstract GaussianProcessKernel factory(String name, List<Parameter> parameters);
        public String toString() { return name; }
        private final String name;
    }

    static GaussianProcessKernel factory(String type, String name, List<Parameter> parameters)
            throws IllegalArgumentException {

        for (AllKernels kernel : AllKernels.values()) {
            if (type.equalsIgnoreCase(kernel.toString())) {
                return kernel.factory(name, parameters);
            }
        }
        throw new IllegalArgumentException("Unknown kernel type");
    }

    abstract class Base extends AbstractModel implements GaussianProcessKernel {

        final List<Parameter> parameters;

        public Base(String name,
                    List<Parameter> parameters) {

            super(name);
            this.parameters = parameters;

            for (Parameter p :parameters) {
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

        @Override
        public double getScale() {
            return parameters.get(0).getParameterValue(0);
        }

        double getLength() {
            return parameters.get(1).getParameterValue(0);
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
