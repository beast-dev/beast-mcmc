package dr.evomodel.substmodel;

import dr.inference.loggers.LogColumn;
import dr.inference.model.Likelihood;
import dr.inference.model.Model;

public interface LogAdditiveCtmcRateProvider extends Model, Likelihood {

    double[] getXBeta();

    LogColumn[] getColumns();

    default double[] getRates() {
        double[] rates = getXBeta();
        for (int i = 0; i < rates.length; ++i) {
            rates[i] = Math.exp(rates[i]);
        }
        return rates;
    }
}
