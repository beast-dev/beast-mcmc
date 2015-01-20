package dr.app.beagle.evomodel.branchmodel.lineagespecific;

import java.util.ArrayList;
import java.util.List;

import dr.app.bss.Utils;
import dr.inference.distribution.ParametricMultivariateDistributionModel;
import dr.inference.model.AbstractModelLikelihood;
import dr.inference.model.CompoundParameter;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;
import dr.inference.model.Variable.ChangeType;

@SuppressWarnings("serial")
public class DirichletProcessPrior  extends AbstractModelLikelihood  {

	private static boolean VERBOSE = false;
	
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
//		this.categoryCount=Utils.findMaximum(categoriesParameter.getParameterValues()) + 1;
        this.N = categoriesParameter.getDimension();
		
		cachedLogFactorials = new ArrayList<Double>();
		cachedLogFactorials.add(0, 0.0);

		// add all
		this.addVariable(this.categoriesParameter);
		this.addVariable(this.gamma);
		
		this.addVariable(this.uniquelyRealizedParameters);
		
		if(baseModel != null) {
		this.addModel(baseModel);
		}
		
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

	public double getGamma() {
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
		return baseModel.logPdf(value);
	}

	public double getRealizedValuesLogDensity() {

		int counts[] = getCounts();
		double total = 0.0;

		//TODO
//		Utils.printArray(counts);
		
		for (int i = 0; i < categoryCount; i++) {

			Parameter param = uniquelyRealizedParameters.getParameter(i);
			int whichCluster = i;

			total += counts[whichCluster] * getLogDensity(param);

			//TODO
//			System.out.println("[" + param.getParameterValue(0) + ", "
//					+ baseModel.getMean()[0] + ", "
//					+ baseModel.getScaleMatrix()[0][0] + "]"
//					+ getLogDensity(param) + " * " + counts[whichCluster]);
			
		}
		
//		System.out.println("total =" + total);
//		System.exit(0);
		
//		return  total;
		return  0;
	}

	public double getCategoriesLogDensity() {

		int[] counts = getCounts();
		
		if (VERBOSE) {
			Utils.printArray(counts);
		}
		
		double loglike = categoryCount * Math.log(getGamma());

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

		this.fireModelChanged();
		likelihoodKnown = false;
		if (!likelihoodKnown) {
			
			logLikelihood = calculateLogLikelihood();
			likelihoodKnown = true;
		}

//		System.out.println(logLikelihood);
		
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
			
			this.fireModelChanged();
			
		} else if (variable == uniquelyRealizedParameters) {

			likelihoodKnown = false;
			this.fireModelChanged();
			
		} else {

			throw new IllegalArgumentException("Unknown parameter");

		}

	}// END: handleVariableChangedEvent

	public void setVerbose() {
		VERBOSE = true;
	}
	
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

		testDirichletProcess(new double[] { 0, 1, 2 }, 3, 1.0, -Math.log(6.0));
		
		testDirichletProcess(new double[] { 0, 0, 1 }, 3, 1.0, -Math.log(6.0));
		
		testDirichletProcess(new double[] { 0, 1, 2, 3, 4 }, 5, 0.5, -6.851184927493743);
		
	}// END: main

	private static void testDirichletProcess(double[] mapping, int categoryCount,double gamma,
			double expectedLogL) {

		Parameter categoriesParameter = new Parameter.Default(mapping);
		Parameter gammaParameter = new Parameter.Default(gamma);

		CompoundParameter dummy = new CompoundParameter("dummy");

		for (int i = 0; i < categoryCount; i++) {
			dummy.addParameter(new Parameter.Default(1.0));
		}

		DirichletProcessPrior dpp = new DirichletProcessPrior(
				categoriesParameter, dummy, null, gammaParameter);
		dpp.setVerbose();
		
		System.out.println("lnL:          " + dpp.getCategoriesLogDensity());
		System.out.println("expected lnL: " + expectedLogL);
	}// END: testDirichletProcess
	
}// END: class
