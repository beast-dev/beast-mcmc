package dr.evomodel.substmodel;

import dr.inference.loggers.LogColumn;
import dr.inference.model.*;

public interface LogAdditiveCtmcRateProvider extends Model, Likelihood {

    double[] getXBeta();

    Parameter getLogRateParameter();

    LogColumn[] getColumns();

    default double[] getRates() {
        double[] rates = getXBeta();
        for (int i = 0; i < rates.length; ++i) {
            rates[i] = Math.exp(rates[i]);
        }
        return rates;
    }

    interface Integrated extends LogAdditiveCtmcRateProvider { }

    interface DataAugmented extends LogAdditiveCtmcRateProvider {

        Parameter getLogRateParameter();

        class Basic extends AbstractModelLikelihood implements DataAugmented {

            private final Parameter logRateParameter;

            public Basic(String name, Parameter logRateParameter) {
                super(name);
                this.logRateParameter = logRateParameter;

                addVariable(logRateParameter);
            }

            public Parameter getLogRateParameter() { return logRateParameter; }

            @Override
            public double[] getXBeta() { // TODO this function should _not_ exponentiate

                final int fieldDim = logRateParameter.getDimension();
                double[] rates = new double[fieldDim];

                for (int i = 0; i < fieldDim; ++i) {
                    rates[i] = Math.exp(logRateParameter.getParameterValue(i));
                }
                return rates;
            }

            @Override
            protected void handleModelChangedEvent(Model model, Object object, int index) { }

            @Override
            protected void handleVariableChangedEvent(Variable variable, int index,
                                                      Parameter.ChangeType type) { }

            @Override
            protected void storeState() { }

            @Override
            protected void restoreState() { }

            @Override
            protected void acceptState() { }

            @Override
            public Model getModel() { return this; }

            @Override
            public double getLogLikelihood() { return 0; }

            @Override
            public void makeDirty() { }
        }
    }
}
