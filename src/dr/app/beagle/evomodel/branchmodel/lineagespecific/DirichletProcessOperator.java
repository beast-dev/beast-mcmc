package dr.app.beagle.evomodel.branchmodel.lineagespecific;

import org.apache.commons.math.MathException;

import dr.inference.distribution.DistributionLikelihood;
import dr.inference.model.CompoundLikelihood;
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
    private int M;
	
	
	private Parameter zParameter;
	private CountableRealizationsParameter countableRealizationsParameter;
	
	private CompoundLikelihood likelihood;
	
	
	
	public DirichletProcessOperator(DirichletProcessPrior dpp, 
			Parameter zParameter, 
			CountableRealizationsParameter countableRealizationsParameter,
			CompoundLikelihood likelihood,
			int M,
			double weight) {

		this.dpp = dpp;
		this.intensity = dpp.getGamma();
		this.uniqueRealizationCount = dpp.getCategoryCount();
		this.realizationCount = zParameter.getDimension();
		
		this.zParameter = zParameter;
		this.countableRealizationsParameter = countableRealizationsParameter;
		this.likelihood = likelihood;
		
		this.M = M;
		setWeight(weight);

	}// END: Constructor

    public Parameter getParameter() {
        return zParameter;
    }//END: getParameter
	
    public Variable getVariable() {
        return zParameter;
    }//END: getVariable
	
	@Override
	public double doOperation() throws OperatorFailedException {

		try {

			doOperate();

		} catch (MathException e) {
			e.printStackTrace();
		}// END: try-catch block

		return 0.0;
	}// END: doOperation

	private void doOperate() throws MathException {

//		int index = 0;
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
			
			double[] clusterProbs = new double[uniqueRealizationCount];
			for (int i = 0; i < uniqueRealizationCount; i++) {

				double logprob = 0;
				if (occupancy[i] == 0) {// draw new

					// draw from base model, evaluate at likelihood
					// M-H for poor people
					
					double loglike = 0.0;
					double candidate = 0.0;
					for (int j = 0; j < M; j++) {
						 candidate = dpp.baseModel.nextRandom()[0];
						 loglike = getPartialLoglike(index, candidate);
					}
					loglike /= M;
					
					logprob = Math.log((intensity) / (realizationCount - 1 + intensity)) + loglike;

				} else {// draw existing
					
					// likelihood for component x_index
					double value = dpp.getUniqueParameter(i).getParameterValue(0);
					double loglike = getPartialLoglike(index, value);
					
					logprob = Math.log(occupancy[i]) / (realizationCount - 1 + intensity) + loglike;

				}// END: occupancy check

				clusterProbs[i] = logprob;
			}// END: i loop

			//rescale
//			double max = dr.app.bss.Utils.max(clusterProbs);
//			for (int i = 0; i < clusterProbs.length; i++) {
//				clusterProbs[i] -=  max;
//			}

			dr.app.bss.Utils.exponentiate(clusterProbs);
//			dr.app.bss.Utils.normalize(clusterProbs);
			
			if (DEBUG) {
				System.out.println("P(z[index] | z[-index]): ");
				dr.app.bss.Utils.printArray(clusterProbs);
			}

			// sample
			int sampledCluster = MathUtils.randomChoicePDF(clusterProbs);
			zParameter.setParameterValue(index, sampledCluster);

			if (DEBUG) {
				System.out.println("sampled category: " + sampledCluster + "\n");
			}

		}// END: index loop

	}// END: doOperate

	private double getPartialLoglike(int index, double candidate) {

		DistributionLikelihood dl = (DistributionLikelihood) likelihood .getLikelihood(index);
		double value = countableRealizationsParameter.getParameterValue(index);

		double loglike = 0.0;
		if (candidate != value) {

			countableRealizationsParameter.setParameterValue(index, candidate);
			loglike = dl.getLogLikelihood();
			countableRealizationsParameter.setParameterValue(index, value);

		} else {

			loglike = dl.getLogLikelihood();

		}

		return loglike;
	}// END: getPartialLoglike
	
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
