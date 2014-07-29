package dr.app.beagle.evomodel.branchmodel.lineagespecific;

import java.util.ArrayList;
import java.util.List;

import dr.app.bss.Utils;
import dr.inference.loggers.LogColumn;
import dr.inference.loggers.Loggable;
import dr.inference.loggers.NumberColumn;
import dr.inference.model.CompoundParameter;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;
import dr.inference.model.Variable.ChangeType;
import dr.inference.model.VariableListener;
import dr.math.distributions.NormalDistribution;

public class DppLogger implements Loggable, VariableListener {

	private Parameter categoriesParameter;
	private CompoundParameter uniquelyRealizedParameters;

	private double[] categoryProbabilities;
	private int newCategoryIndex;
	private double meanForCategory;
	private double newX;

	private int categoryCount;
	private int N;

	public DppLogger(Parameter categoriesParameter, //
			CompoundParameter uniquelyRealizedParameters //
	) {

		this.uniquelyRealizedParameters = uniquelyRealizedParameters;
		this.categoriesParameter = categoriesParameter;
		this.categoryCount = uniquelyRealizedParameters.getDimension();
		this.N = categoriesParameter.getDimension();
		
	}// END: Constructor

	private double[] getCategoryProbs() {

		// int N = categoriesParameter.getDimension();
		double[] categoryProbabilities = new double[categoryCount];

		for (int i = 0; i < N; i++) {
			categoryProbabilities[(int) categoriesParameter
					.getParameterValue(i)]++;
		}// END: N loop

//		Utils.printArray(categoryProbabilities);
//		try {
//		    Thread.sleep(1000);
//		} catch(InterruptedException ex) {
//		    Thread.currentThread().interrupt();
//		}
		
		
		for (int i = 0; i < categoryCount; i++) {
			categoryProbabilities[i] = categoryProbabilities[i] / N;
		}// END: categoryCount loop

		// categoryProbabilities = new double[]{0.10, 0.10, 0.10, 0.10, 0.10,
		// 0.10, 0.10, 0.10, 0.10, 0.10};

		return categoryProbabilities;
	}// END: getCategoryProbs

	@Override
	public LogColumn[] getColumns() {

//		this.categoryProbabilities = getCategoryProbs();

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

		 this.categoryProbabilities = getCategoryProbs();

		this.newCategoryIndex = Utils.sample(categoryProbabilities);
		this.meanForCategory = uniquelyRealizedParameters
				.getParameterValue(newCategoryIndex);

		// baseModel.

		double sd = 1.0;// 0.02857143;
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

	@Override
	public void variableChangedEvent(Variable variable, int index,
			ChangeType type) {
	

		System.out.println("FUBAR");
		
	}

}// END: class
