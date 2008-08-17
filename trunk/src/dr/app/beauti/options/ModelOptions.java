package dr.app.beauti.options;

import dr.app.beauti.priorsPanel.PriorType;
import dr.evolution.util.TaxonList;

import java.util.HashMap;

/**
 * @author Alexei Drummond
 * @author Andrew Rambaut
 */
public class ModelOptions {

	public static final String version = "1.5";
	public static final int YEARS = 0;
	public static final int MONTHS = 1;
	public static final int DAYS = 2;
	public static final int FORWARDS = 0;
	public static final int BACKWARDS = 1;
	public static final int NONE = -1;

	public static final int BIN_SIMPLE = 0;
	public static final int BIN_COVARION = 1;

	public static final int GROWTH_RATE = 0;
	public static final int DOUBLING_TIME = 1;
	public static final int CONSTANT_SKYLINE = 0;
	public static final int LINEAR_SKYLINE = 1;

	public static final int SKYRIDE_UNIFORM_SMOOTHING = 0;
	public static final int SKYRIDE_TIME_AWARE_SMOOTHING = 1;

	public static final int TIME_SCALE = 0;
	public static final int GROWTH_RATE_SCALE = 1;
	public static final int BIRTH_RATE_SCALE = 2;
	public static final int SUBSTITUTION_RATE_SCALE = 3;
	public static final int LOG_STDEV_SCALE = 4;
	public static final int SUBSTITUTION_PARAMETER_SCALE = 5;
	public static final int T50_SCALE = 6;
	public static final int UNITY_SCALE = 7;

	protected void createOperator(String parameterName, OperatorType type, double tuning, double weight) {
		Parameter parameter = getParameter(parameterName);
		operators.put(parameterName, new Operator(parameterName, "", parameter, type, tuning, weight));
	}

	protected void createOperator(String key, String name, String description, String parameterName, OperatorType type, double tuning, double weight) {
		Parameter parameter = getParameter(parameterName);
		operators.put(key, new Operator(name, description, parameter, type, tuning, weight));
	}

	protected void createOperator(String key, String name, String description, String parameterName1, String parameterName2, OperatorType type, double tuning, double weight) {
		Parameter parameter1 = getParameter(parameterName1);
		Parameter parameter2 = getParameter(parameterName2);
		operators.put(key, new Operator(name, description, parameter1, parameter2, type, tuning, weight));
	}

	protected void createScaleOperator(String parameterName, double weight) {
		Parameter parameter = getParameter(parameterName);
		operators.put(parameterName, new Operator(parameterName, "", parameter, OperatorType.SCALE, 0.75, weight));
	}

	protected Parameter createParameter(String name, String description) {
		final Parameter parameter = new Parameter(name, description);
		parameters.put(name, parameter);
		return parameter;
	}

	protected Parameter createParameter(String name, String description, int scale, double value, double lower, double upper) {
		final Parameter parameter = new Parameter(name, description, scale, value, lower, upper);
		parameters.put(name, parameter);
		return parameter;
	}

	protected void createParameter(String name, String description, boolean isNodeHeight, double value, double lower, double upper) {
		parameters.put(name, new Parameter(name, description, isNodeHeight, value, lower, upper));
	}

	protected void createScaleParameter(String name, String description, int scale, double value, double lower, double upper) {
		Parameter p = createParameter(name, description, scale, value, lower, upper);
		p.priorType = PriorType.JEFFREYS_PRIOR;
	}

	protected Parameter createStatistic(String name, String description, boolean isDiscrete) {
		final Parameter parameter = new Parameter(name, description, isDiscrete);
		parameters.put(name, parameter);
		return parameter;
	}

	protected void createStatistic(String name, String description, double lower, double upper) {
		parameters.put(name, new Parameter(name, description, lower, upper));
	}

	public Parameter getParameter(String name) {
		Parameter parameter = parameters.get(name);
		if (parameter == null) {
			for (String key : parameters.keySet()) {
				System.err.println(key);
			}
			throw new IllegalArgumentException("Parameter with name, " + name + ", is unknown");
		}
		return parameter;
	}

	public Parameter getStatistic(TaxonList taxonList) {
		Parameter parameter = statistics.get(taxonList);
		if (parameter == null) {
			for (TaxonList key : statistics.keySet()) {
				System.err.println("Taxon list: " + key.getId());
			}
			throw new IllegalArgumentException("Statistic for taxon list, " + taxonList.getId() + ", is unknown");
		}
		return parameter;
	}

	Operator getOperator(String name) {
		Operator operator = operators.get(name);
		if (operator == null) throw new IllegalArgumentException("Operator with name, " + name + ", is unknown");
		return operator;
	}

	HashMap<String, Parameter> parameters = new HashMap<String, Parameter>();
	HashMap<TaxonList, Parameter> statistics = new HashMap<TaxonList, Parameter>();
	HashMap<String, Operator> operators = new HashMap<String, Operator>();
}
