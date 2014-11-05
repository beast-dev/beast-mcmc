package dr.app.beagle.evomodel.branchmodel.lineagespecific;

import org.apache.commons.math.MathException;

import dr.inference.model.Parameter;
import dr.inference.operators.CoercionMode;
import dr.inference.operators.OperatorFailedException;
import dr.inference.operators.SimpleMCMCOperator;
import dr.math.MathUtils;
import dr.math.distributions.TDistribution;

public class DirichletProcessOperator extends SimpleMCMCOperator
// AbstractCoercableOperator
{

	private static final boolean DEBUG = false;

	private DirichletProcessPrior dpp;

	private Parameter zParameter;
	private int realizationCount;
	private int uniqueRealizationCount;
	private double intensity;

	public DirichletProcessOperator(DirichletProcessPrior dpp,
			Parameter zParameter, int uniqueRealizationCount, double weight) {
		this(dpp, zParameter, null, uniqueRealizationCount, weight);
	}// END: Constructor

	public DirichletProcessOperator(DirichletProcessPrior dpp,
			Parameter zParameter, CoercionMode mode,
			int uniqueRealizationCount, double weight) {

		// super(mode);

		this.dpp = dpp;
		this.intensity = dpp.getGamma();
		this.zParameter = zParameter;
		this.uniqueRealizationCount = uniqueRealizationCount;
		this.realizationCount = zParameter.getDimension();

		setWeight(weight);

	}// END: Constructor

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

		// zParameter.setParameterValue(1, 1);
		// zParameter.setParameterValue(0, 4);
		// zParameter.setParameterValue(2, 4);

//		int index = MathUtils.nextInt(realizationCount);
		
		for (int index = 0; index < realizationCount; index++) {

			// if z[index] is currently a singleton remove
			int zIndex = (int) zParameter.getParameterValue(index);
			int occupied = 0;
			for (int j = 0; j < realizationCount; j++) {

				if (zIndex == zParameter.getParameterValue(j)) {
					occupied++;
				}// END: value check

			}// END: z loop

			if (DEBUG) {
				System.out.println("index: " + index + " value: " + zIndex + " occupancy: " + occupied);
			}

			if (occupied == 1) {

				for (int j = 0; j < zParameter.getDimension(); j++) {

					int zj = (int) zParameter.getParameterValue(j);
					if (zj > zIndex) {
						zParameter.setParameterValue(j, zj - 1);
					}
				}// END: z loop

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

				double prob = 0;
				if (occupancy[i] == 0) {
					// TODO: M-H step here if not conjugate

					double center = 0;
					double scale = 1;
					double df = realizationCount - 1;

					Parameter param = dpp.getUniqueParameter((int) zParameter.getParameterValue(index));
					TDistribution t = new TDistribution(center, scale, df);
					double predDensity = t.pdf(param.getParameterValue(0));

					// draw new
					prob = intensity / (realizationCount - 1 + intensity);
					clusterProbs[i] = Math.log(prob);

				} else {

					Parameter param = dpp.getUniqueParameter((int) zParameter.getParameterValue(index));
					double density = dpp.getLogDensity(param);

					// draw existing
					prob = occupancy[i] / (realizationCount - 1 + intensity);
					clusterProbs[i] = Math.log(prob);

				}

				sum+=prob;
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

//			 printZ();
//			 System.exit(-1);
		}// END: realizations loop


//		 printZ();
	}// END: doOperate

	// @Override
	// public double getCoercableParameter() {
	// // TODO Auto-generated method stub
	// return 0;
	// }
	//
	// @Override
	// public void setCoercableParameter(double value) {
	// // TODO Auto-generated method stub
	//
	// }
	//
	// @Override
	// public double getRawParameter() {
	// // TODO Auto-generated method stub
	// return 0;
	// }

	@Override
	public String getPerformanceSuggestion() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public String getOperatorName() {
		return zParameter.getParameterName();
	}

	private void printZ() {

		for (int i = 0; i < zParameter.getDimension(); i++) {
			System.out.print(zParameter.getParameterValue(i) + " ");
		}
		System.out.println();

	}// END: printZ

}// END: class
