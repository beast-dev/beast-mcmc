package dr.math.distributions.gp;

import dr.inference.model.AbstractModel;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;

import java.util.List;

import static dr.math.ModifiedBesselFirstKind.bessi;
import static dr.math.functionEval.GammaFunction.factorial;
import static dr.math.functionEval.GammaFunction.gamma;

/**
 * @author Marc A. Suchard
 * @author Filippo Monti
 */
public interface GaussianProcessKernel {

    double getCorrelation(double x, double y);

    @SuppressWarnings("unused")
    double getCorrelation(double[] x, double[] y);

    double getScale();

    class DotProduct extends Base {

        public DotProduct(String name, List<Parameter> parameters) {
            super(name, parameters);
        }

        public double getCorrelation(double x, double y) {
            return  x * y;
        }

        public double getCorrelation(double[] x, double[] y) {
            final int dim = x.length;

            double product = 0.0;
            for (int i = 0; i < dim; ++i) {
                product += x[i] * y[i];
            }

            return product;
        }

        private static final String TYPE = "DotProduct";

    }

    class RadialBasisFunction extends Base {

        public RadialBasisFunction(String name, List<Parameter> parameters) {
            super(name, parameters);
        }

        private double functionalForm(double normSquared) {
            double length = parameters.get(0).getParameterValue(0);
            return Math.exp(-normSquared / (2 * length * length));
        }

        public double getCorrelation(double x, double y) {
            double diff = x - y;
            return functionalForm(diff * diff);
        }

        public double getCorrelation(double[] x, double[] y) {
            return functionalForm(normSquared(x, y));
        }

        private static final String TYPE = "RadialBasisFunction";
    }

    class OrnsteinUhlenbeck extends L1Base {

        public OrnsteinUhlenbeck(String name, List<Parameter> parameters) { super(name, parameters); }

        double functionalForm(double norm) {
            double length = parameters.get(0).getParameterValue(0);
            return Math.exp(-norm / length);
        }

        private static final String TYPE = "OrnsteinUhlenbeck";
    }

    class Matern extends L1Base {

        private final Parameter orderParameter;
        private final int order;
        private final double normalization;
        private final double scale;

        public Matern(String name, List<Parameter> parameters) {
            super(name, parameters);
            this.orderParameter = parameters.get(1);

            this.order = (int) orderParameter.getParameterValue(0);
            this.normalization = Math.pow(2, 1 - order) / factorial(order - 1);
            this.scale = Math.sqrt(2 * order);
        }

        double functionalForm(double norm) {
            double length = parameters.get(0).getParameterValue(0);
            double argument = scale * norm / length;

            return normalization * Math.pow(argument, order) * bessi(argument, order);
        }

        @Override
        protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
            if (variable == orderParameter) {
                throw new RuntimeException("Not yet implemented");
            }
            super.handleVariableChangedEvent(variable, index, type);
        }

        private static final String TYPE = "Matern";
    }

    abstract class L1Base extends Base {

        public L1Base(String name, List<Parameter> parameters) { super(name, parameters); }

        abstract double functionalForm(double norm);

        public double getCorrelation(double x, double y) {
            double norm = Math.abs(x - y);
            return functionalForm(norm);
        }

        public double getCorrelation(double[] x, double[] y) {
            double norm = Math.sqrt(normSquared(x, y));
            return functionalForm(norm);
        }
    }

    static GaussianProcessKernel factory(String type, String name, List<Parameter> parameters)
            throws IllegalArgumentException {
        if (type.equalsIgnoreCase(DotProduct.TYPE)) {
            return new DotProduct(name, parameters);
        } else if (type.equalsIgnoreCase(RadialBasisFunction.TYPE)) {
            return new RadialBasisFunction(name, parameters);
        } else if (type.equalsIgnoreCase(Matern.TYPE)) {
            return new Matern(name, parameters);
        } else if (type.equalsIgnoreCase(OrnsteinUhlenbeck.TYPE)) {
            return new OrnsteinUhlenbeck(name, parameters);
        } else {
            throw new IllegalArgumentException("Unknown kernel type");
        }
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
            return 1.0;
            // return parameters.get(0).getParameterValue(0); // TODO
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
