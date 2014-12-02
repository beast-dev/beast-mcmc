package dr.app.beagle.evomodel.branchmodel.lineagespecific;

import org.apache.commons.math.MathException;

import dr.inference.model.Parameter;
import dr.inference.model.Variable;
import dr.inference.operators.GibbsOperator;
import dr.inference.operators.OperatorFailedException;
import dr.inference.operators.SimpleMCMCOperator;
import dr.math.MathUtils;

public class DirichletProcessOperator extends SimpleMCMCOperator implements
		GibbsOperator {

	private static final boolean DEBUG = false;

	private DirichletProcessPrior dpp;

	private Parameter zParameter;
	private int realizationCount;
	private int uniqueRealizationCount;
	private double intensity;

	public DirichletProcessOperator(DirichletProcessPrior dpp, 
			Parameter zParameter, 
			double weight) {

		this.dpp = dpp;
		this.intensity = dpp.getGamma();
		this.zParameter = zParameter;
		this.uniqueRealizationCount = dpp.getCategoryCount();
		this.realizationCount = zParameter.getDimension();

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

		return 0;
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

			if (occupied == 1) {

//				for (int j = 0; j < zParameter.getDimension(); j++) {
//
//					int zj = (int) zParameter.getParameterValue(j);
//					if (zj > zValue) {
//						zParameter.setParameterValue(j, zj - 1);
//					}
//				}// END: z loop

				zParameter.setParameterValue(index, 0);

			}// END: singleton check

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
			double sum = 0;
			for (int i = 0; i < uniqueRealizationCount; i++) {

				double loglike = dpp.getRealizedValuesLogDensity();

				double prob = 0;
				if (occupancy[i] == 0) {// draw new
					
					prob = (intensity / (realizationCount - 1 + intensity));
					
				} else {// draw existing

					prob = (occupancy[i] / (realizationCount - 1 + intensity)) ;

				}//END: occupancy check

//				TODO
//				System.out.println("i: ");
//				System.out.println(prob);
//				System.out.println(Math.log(prob));
//				System.out.println(Math.log(prob) + loglike);
//				System.out.println();
				
				clusterProbs[i] = Math.log(prob) + loglike;
				sum += (prob);
			}// END: i loop

			
			// normalize (b in Neal 2000)
			double logsum = Math.log(sum);
			for (int i = 0; i < clusterProbs.length; i++) {
				clusterProbs[i] -=  logsum;
			}
			
			 dr.app.bss.Utils.exponentiate(clusterProbs);
			
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

	@Override
	public String getOperatorName() {
		return DirichletProcessOperatorParser.DIRICHLET_PROCESS_OPERATOR;
	}

//	private void printZ() {
//		for (int i = 0; i < zParameter.getDimension(); i++) {
//			System.out.print(zParameter.getParameterValue(i) + " ");
//		}
//		System.out.println();
//	}// END: printZ

	@Override
	public String getPerformanceSuggestion() {
		return null;
	}// END: getPerformanceSuggestion

	@Override
	public int getStepCount() {
		return realizationCount;
	}// END: getStepCount

}// END: class
