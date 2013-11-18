package dr.evomodel.epidemiology.casetocase;

import dr.inference.distribution.NormalDistributionModel;
import dr.inference.model.Parameter;
import dr.math.distributions.Distribution;
import dr.math.distributions.NormalDistribution;
import dr.xml.*;

/**
 * Created with IntelliJ IDEA.
 * User: mhall
 * Date: 04/11/2013
 * Time: 14:50
 * To change this template use File | Settings | File Templates.
 */

public abstract class AnalyticallySolvablePosteriorFunction {


    public abstract double getPosteriorValue(Parameter params, double[] data);

    public abstract double getLogPosteriorValue(Parameter params, double[] data);

    public abstract double getPriorValue(Parameter params);

    public abstract double getPriorPredictivePDF(double value);

    public abstract double getPriorPredictiveCDF(double value);

    public abstract double getPriorPredictiveInterval(double start, double end);

    public abstract Distribution getPriorPredictiveDistribution();

}
