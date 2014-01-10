package dr.evomodel.epidemiology.casetocase;

import dr.evolution.coalescent.DemographicFunction;
import dr.evolution.util.Units;
import dr.evomodel.coalescent.DemographicModel;
import dr.inference.model.Parameter;

/**
 * Shifted, four-parameter logistic growth function.
 *
 * To avoid numerical errors, we assume N0>2*N50, and parametrise based on positive values of N50 and N0-2*N50.
 *
 * @author Matthew Hall
 */
public class LogisticGrowthN0N50Model extends DemographicModel {

    public LogisticGrowthN0N50Model(Parameter n50Parameter, Parameter n0Minus2n50Parameter,
                                    Parameter growthRateParameter, Parameter t50Parameter, Units.Type units) {
        this(LogisticGrowthN0N50ModelParser.LOGISTIC_GROWTH_MODEL, n50Parameter, n0Minus2n50Parameter,
                growthRateParameter, t50Parameter, units);
    }

    public LogisticGrowthN0N50Model(String name, Parameter n50Parameter, Parameter n0Minus2n50Parameter,
                                    Parameter growthRateParameter, Parameter t50Parameter, Units.Type units) {
        super(name);

        logisticGrowthNLimit = new LogisticGrowthN0N50(units);

        this.n50Parameter = n50Parameter;
        addVariable(n50Parameter);
        n50Parameter.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, 1));

        this.n0Minus2n50Parameter = n0Minus2n50Parameter;
        addVariable(n0Minus2n50Parameter);
        n50Parameter.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, 1));

        this.growthRateParameter = growthRateParameter;
        addVariable(growthRateParameter);
        growthRateParameter.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, 1));

        this.t50Parameter = t50Parameter;
        addVariable(t50Parameter);
        t50Parameter.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, 1));

        setUnits(units);
    }

    public DemographicFunction getDemographicFunction() {
        logisticGrowthNLimit.setN50(n50Parameter.getParameterValue(0));
        logisticGrowthNLimit.setGrowthRate(growthRateParameter.getParameterValue(0));
        logisticGrowthNLimit.setT50(t50Parameter.getParameterValue(0));
        logisticGrowthNLimit.setN0(n0Minus2n50Parameter.getParameterValue(0) + 2*n50Parameter.getParameterValue(0));

        return logisticGrowthNLimit;
    }

    private LogisticGrowthN0N50 logisticGrowthNLimit;
    private Parameter n50Parameter;
    private Parameter n0Minus2n50Parameter;
    private Parameter growthRateParameter;
    private Parameter t50Parameter;

}
