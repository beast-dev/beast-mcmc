package dr.app.beagle.evomodel.branchmodel.lineagespecific;

import java.util.LinkedHashMap;

import org.apache.commons.math.MathException;

import dr.inference.model.CompoundLikelihood;
import dr.inference.model.Likelihood;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;
import dr.inference.operators.GibbsOperator;
import dr.inference.operators.OperatorFailedException;
import dr.inference.operators.SimpleMCMCOperator;
import dr.math.MathUtils;

@SuppressWarnings("serial")
public class DirichletProcessOperator extends SimpleMCMCOperator implements
		GibbsOperator {

	private static final boolean DEBUG = false;

	private DirichletProcessPrior dpp;

	private int realizationCount;
	private int uniqueRealizationCount;
	private double intensity;
	private int mhSteps;

	private Parameter zParameter;
	// private CountableRealizationsParameter countableRealizationsParameter;
	private Parameter parameter;

	private CompoundLikelihood likelihood;

	public DirichletProcessOperator(DirichletProcessPrior dpp, //
			Parameter zParameter, //
			// CountableRealizationsParameter countableRealizationsParameter,
			Parameter parameter, //
			CompoundLikelihood likelihood, //
			int mhSteps, //
			double weight//
	) {

		this.dpp = dpp;
		this.intensity = dpp.getGamma();
		this.uniqueRealizationCount = dpp.getCategoryCount();
		this.realizationCount = zParameter.getDimension();

		this.zParameter = zParameter;
		// this.countableRealizationsParameter = countableRealizationsParameter;
		this.parameter = parameter;
		this.likelihood = likelihood;

		this.mhSteps = mhSteps;
		setWeight(weight);

	}// END: Constructor

	public Parameter getParameter() {
		return zParameter;
	}// END: getParameter

	public Variable getVariable() {
		return zParameter;
	}// END: getVariable

	@Override
	public double doOperation() throws OperatorFailedException {

		try {

//			doOperate();
			doOp();

		} catch (MathException e) {
			e.printStackTrace();
		}// END: try-catch block

		return 0.0;
	}// END: doOperation

	private void doOp() throws MathException {
		
		for (int index = 0; index < realizationCount; index++) {
		
			int[] occupancy = new int[uniqueRealizationCount];
			for (int i = 0; i < realizationCount; i++) {
				if (i != index) {
					int j = (int) zParameter.getParameterValue(i);
					occupancy[j]++;
				}// END: i check
			}// END: i loop
			
			

			// TODO: set parameters at index values 
			int category = (int) zParameter.getParameterValue(index);
			
			for (int i = 0; i < uniqueRealizationCount; i++) {
			
				double candidate =0;
				if (occupancy[i] == 0) {// draw new

					// draw from base model

					 candidate = dpp.baseModel.nextRandom()[0];

				} else {// draw existing

					// likelihood for component x_index

					 candidate = dpp.getUniqueParameter(i)
							.getParameterValue(0);


				}// END: occupancy check
				
				parameter.setParameterValue(category, candidate);
				
				
			}// END: i loop
			
			double loglike = likelihood.getLogLikelihood();
			System.out.println(loglike);
			System.exit(-1);
			
		}//END: index loop
		
		
		
		
	}//END: doOp
	
	
	private void doOperate() throws MathException {

		// int index = 0;
		for (int index = 0; index < realizationCount; index++) {

			int[] occupancy = new int[uniqueRealizationCount];
			for (int i = 0; i < realizationCount; i++) {

				if (i != index) {

					int j = (int) zParameter.getParameterValue(i);
					occupancy[j]++;

				}// END: i check

			}// END: i loop

			if (DEBUG) {
				System.out.println("N[-index]: ");
				dr.app.bss.Utils.printArray(occupancy);
			}

			Likelihood clusterLikelihood = (Likelihood) likelihood.getLikelihood(index);
//			Likelihood clusterLikelihood = likelihood;
			
			int category = (int) zParameter.getParameterValue(index);
			double value = parameter.getParameterValue(category);
			
			double[] clusterProbs = new double[uniqueRealizationCount];
			
			for (int i = 0; i < uniqueRealizationCount; i++) {

				double logprob = 0;
				if (occupancy[i] == 0) {// draw new

					// draw from base model, evaluate at likelihood

					double candidate = dpp.baseModel.nextRandom()[0];

					parameter.setParameterValue(category, candidate);
					double loglike = clusterLikelihood.getLogLikelihood();
					parameter.setParameterValue(category, value);

					logprob = Math.log((intensity)
							/ (realizationCount - 1 + intensity))
							+ loglike;

				} else {// draw existing

					// likelihood for component x_index

					double candidate = dpp.getUniqueParameter(i)
							.getParameterValue(0);

					parameter.setParameterValue(category, candidate);
					double loglike = clusterLikelihood.getLogLikelihood();
					parameter.setParameterValue(category, value);

					logprob = Math.log(occupancy[i])
							/ (realizationCount - 1 + intensity) + loglike;

				}// END: occupancy check

				clusterProbs[i] = logprob;
			}// END: i loop

			
			//////////////////////////////////////
			
			dr.app.bss.Utils.exponentiate(clusterProbs);

//			dr.app.bss.Utils.printArray(clusterProbs);
//			System.exit(-1);

			if (DEBUG) {
				System.out.println("P(z[index] | z[-index]): ");
				dr.app.bss.Utils.printArray(clusterProbs);
			}

			// sample
			int sampledCluster = MathUtils.randomChoicePDF(clusterProbs);
			zParameter.setParameterValue(index, sampledCluster);

			if (DEBUG) {
				System.out
						.println("sampled category: " + sampledCluster + "\n");
			}

		}// END: index loop

	}// END: doOperate

	@Override
	public String getOperatorName() {
		return DirichletProcessOperatorParser.DIRICHLET_PROCESS_OPERATOR;
	}

	@Override
	public String getPerformanceSuggestion() {
		return null;
	}// END: getPerformanceSuggestion

	@Override
	public int getStepCount() {
		return realizationCount;
	}// END: getStepCount

}// END: class
