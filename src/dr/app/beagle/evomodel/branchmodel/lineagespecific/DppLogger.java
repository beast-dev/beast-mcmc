package dr.app.beagle.evomodel.branchmodel.lineagespecific;

import dr.inference.loggers.LogColumn;
import dr.inference.loggers.Loggable;
import dr.inference.loggers.NumberColumn;
import dr.inference.model.CompoundParameter;
import dr.inference.model.Parameter;
import dr.math.distributions.NormalDistribution;

public class DppLogger implements Loggable {

	public enum LogMode {
		NEW_VALUE, NEW_CATEGORY, CATEGORY_PROBABILITIES;
	}

	private Parameter categoryProbabilitiesParameter;
	private CompoundParameter uniquelyRealizedParameters;
	private LogMode mode;

	private double[] probs;
	private int categoryIndex;
	private double meanForCategory;
	private double newX;

	public DppLogger(Parameter categoryProbabilitiesParameter,
			CompoundParameter uniquelyRealizedParameters, LogMode mode) {

		this.categoryProbabilitiesParameter = categoryProbabilitiesParameter;
		this.uniquelyRealizedParameters = uniquelyRealizedParameters;
		this.mode = mode;

	}// END: Constructor

	@Override
	public LogColumn[] getColumns() {

		doStuff();

		LogColumn[] columns = null;
		if (mode == LogMode.NEW_VALUE) {

			columns = new LogColumn[] { new NewLogger() };

		} else if (mode == LogMode.NEW_CATEGORY) {

			columns = new LogColumn[] { new NewCategoryLogger() };

		} else if (mode == LogMode.CATEGORY_PROBABILITIES) {

			int dimension = categoryProbabilitiesParameter.getDimension();
			columns = new LogColumn[dimension];
			for (int i = 0; i < dimension; i++) {

				columns[i] = new ProbabilitiesLogger(i);

			}

		} else {
			// do nothing
		}

		return columns;
	}// END: getColumns

	private void doStuff() {

		this.probs = categoryProbabilitiesParameter.getParameterValues();
		this.categoryIndex = dr.app.bss.Utils.sample(probs);
		this.meanForCategory = uniquelyRealizedParameters
				.getParameterValue(categoryIndex);

		double sd = 1.0;
		NormalDistribution nd = new NormalDistribution(meanForCategory, sd);
		this.newX = (double) nd.nextRandom();

	}

	private class NewLogger extends NumberColumn {

		public NewLogger() {
			super("x.new");
		}

		@Override
		public double getDoubleValue() {

			doStuff();

			return newX;
		}// END: getDoubleValue

	}// END: NewLogger class

	private class ProbabilitiesLogger extends NumberColumn {

		private int i;

		public ProbabilitiesLogger(int i) {
			super("category.prob." + i);
			this.i = i;
		}

		@Override
		public double getDoubleValue() {

			doStuff();

			return probs[i];
		}// END: getDoubleValue

	}// END: NewCategoryLogger class

	private class NewCategoryLogger extends NumberColumn {

		public NewCategoryLogger() {
			super("category.new");
		}

		@Override
		public double getDoubleValue() {

			doStuff();

			return categoryIndex;
		}// END: getDoubleValue

	}// END: NewCategoryLogger class

}// END: class
