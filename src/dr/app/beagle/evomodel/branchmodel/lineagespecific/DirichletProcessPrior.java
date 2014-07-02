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
public class DirichletProcessPrior extends AbstractModelLikelihood {

	private Parameter categoriesParameter;
	private ArrayList<Integer> mapping;
	private int mappingLength;
	private int categoryCount;
	private boolean likelihoodKnown = false;
	private double logLikelihood;
	private boolean mappingsKnown = false;
	private CompoundParameter uniquelyRealizedParameters;
	private ParametricMultivariateDistributionModel baseModel;
	private Parameter gamma;
	private int[] counts;
	private final List<Double> cachedLogFactorials;

	public DirichletProcessPrior(Parameter categoriesParameter, //
			CompoundParameter uniquelyRealizedParameters, //
			ParametricMultivariateDistributionModel baseModel, //
			Parameter gamma //
	) {

		super("");

		this.categoriesParameter = categoriesParameter;
		getMappings();

		this.baseModel = baseModel;

		this.uniquelyRealizedParameters = uniquelyRealizedParameters;

		// M clusters
		this.categoryCount = uniquelyRealizedParameters.getDimension();
		// vector z of cluster assignments
		this.mappingLength = mapping.size();

		this.gamma = gamma;

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

		//TODO
//		System.out.println("i:" + i);
//		System.out.println("size:" + cachedLogFactorials.size());
		
		if ( cachedLogFactorials.size() <= i) {

//			//TODO
//			System.out.println("FUBAR");
			
			for (int j = cachedLogFactorials.size() - 1; j <= i; j++) {

				double logfactorial = cachedLogFactorials.get(j) + Math.log(j + 1);
				cachedLogFactorials.add(logfactorial);
			}

		}

//		System.out.println(i);
		
		return cachedLogFactorials.get(i);
	}

	private void getMappings() {

		if (!mappingsKnown) {

			mapping = new ArrayList<Integer>();
			for (int i = 0; i < categoriesParameter.getDimension(); i++) {

				int category = (int) categoriesParameter.getParameterValue(i);
				mapping.add(i, category);

			}//END: i loop
			
			mappingsKnown = true;
		}// END: mappingsKnown check
		
	}// END: getMappings

	/**
	 * Assumes mappings start from index 0
	 * */
	private int[] getCounts() {

		// eta_k parameters (number of assignments to each category)
		counts = new int[categoryCount];
		for (int i = 0; i < mappingLength; i++) {

			counts[mapping.get(i)]++;

		}// END: i loop

		return counts;
	}// END: getCounts

	private int getMapping(int i) {
		return mapping.get(i);
	}

	private int getCount(int whichDensity) {
		return counts[whichDensity];
	}

	private double getLogDensity(Parameter parameter) {
		return baseModel.logPdf(parameter.getAttributeValue());
		
	}

	private double getRealizedValuesLogDensity() {

		double total = 0.0;

		for (int i = 0; i < uniquelyRealizedParameters.getParameterCount(); ++i) {

			Parameter param = uniquelyRealizedParameters.getParameter(i);
			int whichDensity = getMapping(i);

//			System.out.println("whichDensity:" + whichDensity);
//			System.exit(-1);
			
			
			total += getCount(whichDensity) * getLogDensity(param);

		}

		return total;
	}

	private double getCategoriesLogDensity() {

		int[] counts = getCounts();
		double loglike = categoryCount * Math.log(gamma.getParameterValue(0));

		//TODO
//		Utils.printArray(counts);
		
		for (int k = 0; k < categoryCount; k++) {

			int eta = counts[k];

			//TODO
//			System.out.println("eta:" + eta);
			
			if (eta > 0) {
				loglike += getLogFactorial(eta - 1);
			}

		}// END: k loop

		for (int i = 1; i <= mappingLength; i++) {
			loglike -= Math.log(gamma.getParameterValue(0) + i - 1);
		}// END: i loop

		return loglike;
	}// END: getPriorLoglike

	@Override
	public Model getModel() {
		return this;
	}

	@Override
	public double getLogLikelihood() {

//		System.out.println("FUBAR");
//		System.exit(-1);
		
		if (!likelihoodKnown) {
			
			logLikelihood = calculateLogLikelihood();
			likelihoodKnown = true;
		}

		return logLikelihood;
	}

	private double calculateLogLikelihood() {
		
		
		return getCategoriesLogDensity() + getRealizedValuesLogDensity();
	}

	@Override
	public void makeDirty() {
		likelihoodKnown = false;
		mappingsKnown = false;
	}

	@Override
	protected void handleModelChangedEvent(Model model, Object object, int index) {
		likelihoodKnown = false;
	}

	@Override
	protected void handleVariableChangedEvent(Variable variable, int index,
			ChangeType type) {

		if (variable == categoriesParameter) {
			
			mappingsKnown = false;
			this.fireModelChanged();
			
		} else if (variable == gamma) {

			// do nothing

		} else if (variable == uniquelyRealizedParameters) {

			likelihoodKnown = false;
			this.fireModelChanged();
			
			// TODO: do sth

		} else {

			// assuming list<list> hyperparameterList changed

			throw new IllegalArgumentException("Unknown parameter");

		}

		likelihoodKnown = true;

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
		// Parameter[] params = new Parameter[dim*dim];

		MathUtils.setSeed(666);

		for (int i = 0; i < dim; i++) {
			mean[i] = MathUtils.nextDouble();
			for (int j = i; j < dim; j++) {
				precision[i + j * dim] = MathUtils.nextDouble();
			}
		}

		Parameter meanParameter = new Parameter.Default(mean);

		Parameter prec1 = new Parameter.Default(new double[] { precision[0] });
//		Parameter prec2 = new Parameter.Default(new double[] { precision[2],
//				precision[3] });
		Parameter[] params = new Parameter[] { prec1 };
		MatrixParameter precParameter = new MatrixParameter("name", params);

		MultivariateNormalDistributionModel baseModel = new MultivariateNormalDistributionModel(
				meanParameter, precParameter);

		Parameter gamma = new Parameter.Default(0.9);

		DirichletProcessPrior dpp = new DirichletProcessPrior(
				categoriesParameter, realizedThetas, baseModel, gamma);
		System.out.println(dpp.getLogLikelihood());

	}// END: main

}// END: class

