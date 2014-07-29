package dr.app.beagle.evomodel.branchmodel.lineagespecific;

import java.util.ArrayList;
import java.util.List;

import dr.app.bss.Utils;
import dr.inference.distribution.MultivariateNormalDistributionModel;
import dr.inference.distribution.ParametricMultivariateDistributionModel;
import dr.inference.model.AbstractModelLikelihood;
import dr.inference.model.CompoundParameter;
import dr.inference.model.MatrixParameter;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;
import dr.inference.model.Variable.ChangeType;
import dr.math.MathUtils;

@SuppressWarnings("serial")
public class DirichletProcessPrior  extends AbstractModelLikelihood  {

	private Parameter categoriesParameter;
	private CompoundParameter uniquelyRealizedParameters;
	private ParametricMultivariateDistributionModel baseModel;
	private Parameter gamma;

	private int categoryCount;
	private int N;
	
	private boolean likelihoodKnown = false;
	private double logLikelihood;
	private final List<Double> cachedLogFactorials;

	public DirichletProcessPrior(Parameter categoriesParameter, //
			CompoundParameter uniquelyRealizedParameters, //
			ParametricMultivariateDistributionModel baseModel, //
			Parameter gamma //
	) {

		super("");

		// vector z of cluster assignments
		this.categoriesParameter = categoriesParameter;
		this.baseModel = baseModel;
		this.uniquelyRealizedParameters = uniquelyRealizedParameters;
		this.gamma = gamma;
		
		// K clusters
		this.categoryCount = uniquelyRealizedParameters.getDimension();
        this.N = categoriesParameter.getDimension();
		
		cachedLogFactorials = new ArrayList<Double>();
		cachedLogFactorials.add(0, 0.0);

		// add all
		this.addVariable(this.categoriesParameter);
		this.addVariable(this.gamma);
		this.addVariable(this.uniquelyRealizedParameters);

		this.addModel(baseModel);

		this.likelihoodKnown = false;
	}// END: Constructor

	private double getLogFactorial(int i) {

		if ( cachedLogFactorials.size() <= i) {

			for (int j = cachedLogFactorials.size() - 1; j <= i; j++) {

				double logfactorial = cachedLogFactorials.get(j) + Math.log(j + 1);
				cachedLogFactorials.add(logfactorial);
			}

		}

		return cachedLogFactorials.get(i);
	}

	/**
	 * Assumes mappings start from index 0
	 * */
	private int[] getCounts() {

		// eta_k parameters (number of assignments to each category)
		int[] counts = new int[categoryCount];
		for (int i = 0; i < N; i++) {

			int category = getMapping(i);
			counts[category]++;

		}// END: i loop

		return counts;
	}// END: getCounts

	private double getGamma() {
		return gamma.getParameterValue(0);
	}
	
	private int getMapping(int i) {
		return (int) categoriesParameter.getParameterValue(i);
	}

//	private int getCount(int whichCluster) {
//		return counts[whichCluster];
//	}

	private double getLogDensity(Parameter parameter) {
		
		double value[] = parameter.getAttributeValue();
		
//		try {
//			Utils.printArray(value);
//			Thread.sleep(1000);
//		} catch (InterruptedException e) {
//			e.printStackTrace();
//		}
		
		return baseModel.logPdf(value);
		
	}

	private double getRealizedValuesLogDensity() {

		int counts[] = getCounts();
		double total = 0.0;

		for (int i = 0; i < categoryCount; ++i) {

			Parameter param = uniquelyRealizedParameters.getParameter(i);
			int whichCluster = getMapping(i);

			total += counts[whichCluster] * getLogDensity(param);

		}

		return total;
	}

	private double getCategoriesLogDensity() {

		int[] counts = getCounts();
		double loglike = categoryCount * Math.log(gamma.getParameterValue(0));

		for (int k = 0; k < categoryCount; k++) {

			int eta = counts[k];

			if (eta > 0) {
				loglike += getLogFactorial(eta - 1);
			}

		}// END: k loop

		for (int i = 1; i <= N; i++) {
			loglike -= Math.log(getGamma() + i - 1);
		}// END: i loop

		return loglike;
	}// END: getPriorLoglike

	@Override
	public Model getModel() {
		return this;
	}

	@Override
	public double getLogLikelihood() {

		likelihoodKnown = false;
		if (!likelihoodKnown) {
			
			logLikelihood = calculateLogLikelihood();
			likelihoodKnown = true;
		}

		return logLikelihood;
	}

	private double calculateLogLikelihood() {
		return getCategoriesLogDensity() + getRealizedValuesLogDensity();
	}//END: calculateLogLikelihood

	@Override
	public void makeDirty() {
//		likelihoodKnown = false;
	}

	@Override
	protected void handleModelChangedEvent(Model model, Object object, int index) {
		likelihoodKnown = false;
	}

	@Override
	protected void handleVariableChangedEvent(Variable variable, int index,
			ChangeType type) {

		if (variable == categoriesParameter) {
			
			this.fireModelChanged();
			
		} else if (variable == gamma) {

			// do nothing

		} else if (variable == uniquelyRealizedParameters) {

			likelihoodKnown = false;
			this.fireModelChanged();
			
		} else {

			throw new IllegalArgumentException("Unknown parameter");

		}

	}// END: handleVariableChangedEvent

	@Override
	protected void storeState() {

	}

	@Override
	protected void restoreState() {
		likelihoodKnown = false;
	}

	@Override
	protected void acceptState() {

	}

	public static void main(String args[]) {

		double[] mappings = new double[] { 0, 0, 1, 1, 2 };
		Parameter categoriesParameter = new Parameter.Default(mappings);

		Parameter rTheta1 = new Parameter.Default(1.0);
		Parameter rTheta2 = new Parameter.Default(1.0);
		Parameter rTheta3 = new Parameter.Default(2.0);
		Parameter rTheta4 = new Parameter.Default(2.0);
		Parameter rTheta5 = new Parameter.Default(3.0);

		CompoundParameter realizedThetas = new CompoundParameter(
				"realizedThetas");
		realizedThetas.addParameter(rTheta1);
		realizedThetas.addParameter(rTheta2);
		realizedThetas.addParameter(rTheta3);
		realizedThetas.addParameter(rTheta4);
		realizedThetas.addParameter(rTheta5);

		int dim = 1;
		double[] mean = new double[dim];
		double[] precision = new double[dim * dim];

		MathUtils.setSeed(666);

		for (int i = 0; i < dim; i++) {
			mean[i] = MathUtils.nextDouble();
			for (int j = i; j < dim; j++) {
				precision[i + j * dim] = MathUtils.nextDouble();
			}
		}

		Parameter meanParameter = new Parameter.Default(mean);

		Parameter prec1 = new Parameter.Default(new double[] { precision[0] });
		Parameter[] params = new Parameter[] { prec1 };
		MatrixParameter precParameter = new MatrixParameter("precision", params);

		MultivariateNormalDistributionModel baseModel = new MultivariateNormalDistributionModel(
				meanParameter, precParameter);

		Parameter gamma = new Parameter.Default(0.9);

		DirichletProcessPrior dpp = new DirichletProcessPrior(
				categoriesParameter, realizedThetas, baseModel, gamma);
		
		System.out.println("lnL:" + dpp.getLogLikelihood());

	}// END: main

}// END: class

