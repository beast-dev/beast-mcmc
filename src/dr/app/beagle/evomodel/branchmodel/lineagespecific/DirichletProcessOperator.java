package dr.app.beagle.evomodel.branchmodel.lineagespecific;

import org.apache.commons.math.MathException;

import dr.inference.distribution.DistributionLikelihood;
import dr.inference.distribution.ParametricDistributionModel;
import dr.inference.model.CompoundLikelihood;
import dr.inference.model.Likelihood;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;
import dr.inference.operators.GibbsOperator;
import dr.inference.operators.OperatorFailedException;
import dr.inference.operators.SimpleMCMCOperator;
import dr.inferencexml.distribution.DistributionModelParser;
import dr.math.MathUtils;
import dr.math.distributions.NormalDistribution;

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

//		int index = MathUtils.nextInt(realizationCount);
		for (int index = 0; index < realizationCount; index++) {

			// if z[index] is currently a singleton remove
			int zValue = (int) zParameter.getParameterValue(index);
			
			int occupied = 0;
			for (int j = 0; j < realizationCount; j++) {

				if (zValue == zParameter.getParameterValue(j)) {
					occupied++;
				}// END: value check

			}// END: z loop

			if (DEBUG) {
				System.out.println("index: " + index + " value: " + zValue + " occupancy: " + occupied);
			}

//			if (occupied == 1) {
//				for (int j = 0; j < zParameter.getDimension(); j++) {
//					int zj = (int) zParameter.getParameterValue(j);
//					if (zj > zValue) {
//						zParameter.setParameterValue(j, zj - 1);
//					}
//				}// END: z loop
//				zParameter.setParameterValue(index, 0);
//			}// END: singleton check

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

				double prob = 0;
				if (occupancy[i] == 0) {// draw new

					// draw from base model, evaluate at likelihood

					double[] candidate = dpp.baseModel.nextRandom();
					DistributionLikelihood dl = (DistributionLikelihood) likelihood .getLikelihood(index);
					double[] data = dl.getDataList().get(0) .getAttributeValue();
					ParametricDistributionModel dm = (ParametricDistributionModel) likelihood .getModel().getModel(index);

					// TODO: general implementation for every distribution model
					// TODO: move to method
					int m = 10;
					double loglike = 0;
					for (int j = 0; j < m; j++) {

						loglike += NormalDistribution.logPdf(data[0],
								candidate[0], (double) dm.getVariable(1)
										.getValue(0));

					}
					loglike /= m;

					prob = Math.log((intensity)
							/ (realizationCount - 1 + intensity))
							+ loglike;

				} else {// draw existing

					// likelihood for component x_index
					double loglike = likelihood.getLikelihood(index) .getLogLikelihood();
					prob = Math.log((occupancy[i]) / (realizationCount - 1 + intensity)) + loglike;

				}// END: occupancy check

				clusterProbs[i] = prob;
			}// END: i loop

			//rescale
			double max = dr.app.bss.Utils.max(clusterProbs);
			for (int i = 0; i < clusterProbs.length; i++) {
				clusterProbs[i] -=  max;
			}
			
			dr.app.bss.Utils.exponentiate(clusterProbs);
			dr.app.bss.Utils.normalize(clusterProbs);
			
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

//			printZ();
//			System.exit(-1);	
			
	}// END: doOperate

//	private double getPartialLikelihood(int index) {
//		double like = 0.0;
//		
//	Model model = 	likelihood.getModel().getModel(index);
//		
//	
//	return like;
//	}
	
	@Override
	public String getOperatorName() {
		return DirichletProcessOperatorParser.DIRICHLET_PROCESS_OPERATOR;
	}

	private void printZ() {
		for (int i = 0; i < zParameter.getDimension(); i++) {
			System.out.print(zParameter.getParameterValue(i) + " ");
		}
		System.out.println();
	}// END: printZ

	@Override
	public String getPerformanceSuggestion() {
		return null;
	}// END: getPerformanceSuggestion

	@Override
	public int getStepCount() {
		return realizationCount;
	}// END: getStepCount

	private class Dummy {
		  public ParametricDistributionModel dm;

		  public Dummy(ParametricDistributionModel dm) {
		    this.dm = dm; // you can access  
		    
		    
		  }

		  
		  
		  
		}//END: dummy
	
}// END: class
