package dr.app.beagle.evomodel.branchmodel.lineagespecific;

import java.util.ArrayList;
import java.util.List;

import dr.inference.loggers.LogColumn;
import dr.inference.loggers.Loggable;
import dr.inference.loggers.NumberColumn;
import dr.inference.model.CompoundParameter;
import dr.inference.model.Parameter;
import dr.math.distributions.NormalDistribution;

public class DppLogger implements Loggable {

	private Parameter categoryProbabilitiesParameter;
	private CompoundParameter uniquelyRealizedParameters;

	private double[] categoryProbabilities;
	private int newCategoryIndex;
	private double meanForCategory;
	private double newX;

	public DppLogger(Parameter categoryProbabilitiesParameter, //
			CompoundParameter uniquelyRealizedParameters //
	) {

		this.categoryProbabilitiesParameter = categoryProbabilitiesParameter;
		this.uniquelyRealizedParameters = uniquelyRealizedParameters;

	}// END: Constructor

	@Override
	public LogColumn[] getColumns() {

		List<LogColumn> columns = new ArrayList<LogColumn>();
		
		columns.add(new NewLogger("x.new"));
		columns.add(new NewCategoryLogger("category.new"));
		columns.add(new NewMeanLogger("mean.new"));
		for (int i = 0; i < uniquelyRealizedParameters.getDimension(); i++) {
			columns.add(new ProbabilitiesLogger("pi.", i));
		}

		LogColumn[] rtnColumns = new LogColumn[columns.size()];
		
		return columns.toArray(rtnColumns);
	}// END: getColumns

	private void getNew() {

		this.categoryProbabilities = categoryProbabilitiesParameter.getParameterValues();
//		this.categoryProbabilities = new double[]{0.20, 0.20, 0.20, 0.20, 0.20};
		
		this.newCategoryIndex = dr.app.bss.Utils.sample(categoryProbabilities);
		this.meanForCategory = uniquelyRealizedParameters.getParameterValue(newCategoryIndex);

		double sd = 0.02857143;
		NormalDistribution nd = new NormalDistribution(meanForCategory, sd);
		this.newX = (Double) nd.nextRandom();

	}

	private class NewLogger extends NumberColumn {

		public NewLogger(String label) {
			super(label);
		}

		@Override
		public double getDoubleValue() {

			getNew();

			return newX;
		}// END: getDoubleValue

	}// END: NewLogger class

	private class NewMeanLogger extends NumberColumn {

		public NewMeanLogger(String label) {
			super(label);
		}

		@Override
		public double getDoubleValue() {

			return meanForCategory;
		}// END: getDoubleValue

	}// END: NewCategoryLogger class

	private class ProbabilitiesLogger extends NumberColumn {

		private int i;

		public ProbabilitiesLogger(String label, int i) {
			super(label + i);
			this.i = i;
		}

		@Override
		public double getDoubleValue() {

			return categoryProbabilities[i];
		}// END: getDoubleValue

	}// END: NewCategoryLogger class

	private class NewCategoryLogger extends NumberColumn {

		public NewCategoryLogger(String label) {
			super(label);
		}

		@Override
		public double getDoubleValue() {

			return newCategoryIndex;
		}// END: getDoubleValue

	}// END: NewCategoryLogger class

}// END: class
