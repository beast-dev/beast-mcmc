package dr.app.beagle.evomodel.branchmodel.lineagespecific;

import org.apache.commons.math.MathException;

import dr.inference.distribution.DistributionLikelihood;
import dr.inference.distribution.ParametricDistributionModel;
import dr.inference.model.CompoundLikelihood;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;
import dr.inference.operators.GibbsOperator;
import dr.inference.operators.OperatorFailedException;
import dr.inference.operators.SimpleMCMCOperator;
import dr.math.MathUtils;
import dr.math.distributions.NormalDistribution;

@SuppressWarnings("serial")
public class DirichletProcessOperator extends SimpleMCMCOperator implements
		GibbsOperator {

	private static final boolean DEBUG = false;

	private DirichletProcessPrior dpp;

	private Parameter zParameter;
	private int realizationCount;
	private int uniqueRealizationCount;
	private double intensity;

	private CompoundLikelihood likelihood;
	
	public DirichletProcessOperator(DirichletProcessPrior dpp, 
			Parameter zParameter, 
			CompoundLikelihood likelihood,
			double weight) {

		this.dpp = dpp;
		this.intensity = dpp.getGamma();
		this.zParameter = zParameter;
		this.uniqueRealizationCount = dpp.getCategoryCount();
		this.realizationCount = zParameter.getDimension();

		
		this.likelihood = likelihood;
		
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
					
					int m = 1;
					double loglike = 0.0;
					double candidate = 0.0;
					for (int j = 0; j < m; j++) {
						 candidate = dpp.baseModel.nextRandom()[0];
						loglike = getPartialLoglike(index, candidate);
					}
					loglike /= m;
					
					logprob = Math.log((intensity) / (realizationCount - 1 + intensity)) + loglike;

				} else {// draw existing
					
					// likelihood for component x_index
					double mu = dpp.getUniqueParameter(i).getParameterValue(0);
					double loglike = getPartialLoglike(index, mu);
					
					double prob = (occupancy[i]) / (realizationCount - 1 + intensity);
					logprob = Math.log(prob) + loglike;

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

			//TODO
//			System.exit(-1);	
			
		}// END: index loop

	}// END: doOperate

	
	//TODO: general for arbitrary distributions
	private double getPartialLoglike(int index, double parameterValue) {

		DistributionLikelihood dl = (DistributionLikelihood) likelihood .getLikelihood(index);
		ParametricDistributionModel dm = (ParametricDistributionModel) likelihood .getModel().getModel(index);
		double stdev = (Double) dm.getVariable(1) .getValue(0);
		double data = dl.getDataList().get(0) .getAttributeValue()[0];
		
		double loglike = NormalDistribution.logPdf(data, parameterValue, stdev);
		
		return loglike;
	}
	
	
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
