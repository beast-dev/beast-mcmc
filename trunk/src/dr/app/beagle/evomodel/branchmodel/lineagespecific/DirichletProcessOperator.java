package dr.app.beagle.evomodel.branchmodel.lineagespecific;

import java.util.Arrays;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Set;

import org.apache.commons.math.MathException;
import org.apache.commons.math.distribution.TDistributionImpl;

import dr.app.bss.Utils;
import dr.inference.distribution.TDistributionModel;
import dr.inference.model.Parameter;
import dr.inference.operators.AbstractCoercableOperator;
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
			for (int i = 0; i < uniqueRealizationCount; i++) {

				if (occupancy[i] == 0) {
					// TODO: M-H step here if not conjugate

					double center = 0;
					double scale = 1;
					double df = realizationCount - 1;

					Parameter param = dpp.getUniqueParameter((int) zParameter.getParameterValue(index));
					TDistribution t = new TDistribution(center, scale, df);
					double predProb = t.pdf(param.getParameterValue(0));

					// draw new
					clusterProbs[i] = (intensity / (realizationCount - 1 + intensity));

				} else {

					Parameter param = dpp.getUniqueParameter((int) zParameter.getParameterValue(index));
					double prob = dpp.getLogDensity(param);

					// draw existing
					clusterProbs[i] = (occupancy[i] / (realizationCount - 1 + intensity));

				}

			}// END: i loop

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

			// printZ();

		}// END: realizations loop

		// System.exit(-1);

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
