package dr.app.beagle.evomodel.branchmodel.lineagespecific;

import dr.inference.loggers.LogColumn;
import dr.inference.loggers.Loggable;
import dr.inference.loggers.NumberColumn;
import dr.inference.model.CompoundParameter;
import dr.inference.model.Parameter;
import dr.math.distributions.NormalDistribution;

public class DppLogger implements Loggable {

	private Parameter categoryProbabilitiesParameter;
	private CompoundParameter uniquelyRealizedParameters;

	public DppLogger(Parameter categoryProbabilitiesParameter,
			CompoundParameter uniquelyRealizedParameters) {

		this.categoryProbabilitiesParameter = categoryProbabilitiesParameter;
		this.uniquelyRealizedParameters = uniquelyRealizedParameters;

	}// END: Constructor

	@Override
	public LogColumn[] getColumns() {

		int N = 1;
		NewLogger[] columns = new NewLogger[N];

		for (int i = 0; i < N; i++) {

			columns[i] = new NewLogger();

		}

		return columns;
	}// END: getColumns

	private class NewLogger extends NumberColumn {

		public NewLogger() {
			super("x.new");
		}

		@Override
		public double getDoubleValue() {

			double sd = 1.0;

			double[] probs = categoryProbabilitiesParameter
					.getParameterValues();

			int categoryIndex = dr.app.bss.Utils.sample(probs);
			double meanForCategory = uniquelyRealizedParameters
					.getParameterValue(categoryIndex);

//			Utils.printArray(probs);
//			System.out.println("idx: " + categoryIndex + " meanForCat: " + meanForCategory);

			NormalDistribution nd = new NormalDistribution(meanForCategory, sd);

			return (double) nd.nextRandom();
		}// END: getDoubleValue

	}// END: DPColumn class

}// END: class
