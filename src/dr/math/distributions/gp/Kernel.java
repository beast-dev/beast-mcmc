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
public interface Kernel {

    double getCorrelation(double x, double y);

    class DotProduct extends Base implements Kernel {

        public DotProduct(String name, List<Parameter> parameters) {
            super(name, parameters);
        }

        public double getCorrelation(double x, double y) {
            double sigma = parameters.get(0).getParameterValue(0);
            return sigma * x * y;
        }
    }

    class RadialBasisFunction extends Base {

        public RadialBasisFunction(String name, List<Parameter> parameters) {
            super(name, parameters);
        }

        public double getCorrelation(double x, double y) {
            double sigma = parameters.get(0).getParameterValue(0);
            double length = parameters.get(1).getParameterValue(0);
            double diff = x - y;

            return sigma * Math.exp(-(diff * diff) / (2 * length * length));
        }
    }

    class Base extends AbstractModel {

        final List<Parameter> parameters;

        public Base(String name,
                    List<Parameter> parameters) {

            super(name);
            this.parameters = parameters;

            for (Parameter p :parameters) {
                addVariable(p);
            }
        }

        @Override
        protected void handleModelChangedEvent(Model model, Object object, int index) { }

        @Override
        protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
            if (parameters.contains((Parameter) variable)) {
                fireModelChanged(variable, index);
            } else {
                throw new RuntimeException("Unknown variable");
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
