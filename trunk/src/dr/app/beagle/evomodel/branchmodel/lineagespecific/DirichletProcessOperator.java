package dr.app.beagle.evomodel.branchmodel.lineagespecific;

import java.util.ArrayList;

import org.apache.commons.math.MathException;

import dr.inference.model.Parameter;
import dr.inference.operators.AbstractCoercableOperator;
import dr.inference.operators.CoercionMode;
import dr.inference.operators.OperatorFailedException;
import dr.math.MathUtils;
import dr.math.distributions.BetaDistribution;

public class DirichletProcessOperator extends AbstractCoercableOperator {

	private Parameter zParameter;
	private Parameter categoryProbabilitiesParameter;
	private int realizationCount;
	private int uniqueRealizationCount;
	private double intensity;
	private BetaDistribution beta;

	public DirichletProcessOperator(Parameter zParameter, double intensity, int uniqueRealizationCount) {
		this(zParameter, null, null, intensity, uniqueRealizationCount);
	}// END: Constructor

	public DirichletProcessOperator(Parameter zParameter, Parameter categoryProbabilitiesParameter, double intensity, int uniqueRealizationCount) {
		this(zParameter, categoryProbabilitiesParameter, null, intensity, uniqueRealizationCount);
	}// END: Constructor

	public DirichletProcessOperator(Parameter zParameter, Parameter categoryProbabilitiesParameter, CoercionMode mode, double intensity, 
			int uniqueRealizationCount) {

		super(mode);

		this.zParameter = zParameter;
		this.categoryProbabilitiesParameter = categoryProbabilitiesParameter;		
		this.intensity = intensity;
		this.uniqueRealizationCount = uniqueRealizationCount;
		
		realizationCount = zParameter.getDimension();
		beta = new BetaDistribution(this.intensity, 1.0);

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

		beta = new BetaDistribution(0.9, 1.0);
		double[] probs = stickBreak(uniqueRealizationCount);

		if (categoryProbabilitiesParameter != null) {

			for (int i = 0; i < uniqueRealizationCount; i++) {
				double prob = probs[i];
				categoryProbabilitiesParameter.setParameterValue(i, prob);
				
//				System.out.println(prob);
				
			}

		}// END: null check

		for (int i = 0; i < realizationCount; i++) {

			int categoryIndex = dr.app.bss.Utils.sample(probs);
			zParameter.setParameterValue(i, categoryIndex);

		}// END: realizationCount loop

	}// END: doOperate

	public double[] stickBreak(int K) throws MathException {

		double[] probabilities = new double[K];

		double stickLength = 1.0;
		int i = 0;
		for (int x = 1; x <= K; x++) {

			double beta = rBeta();
			probabilities[i] = beta * stickLength;
			stickLength *= (1 - beta);

			i = i + 1;

		}// END: K loop

		double fallThrough = 1 - dr.app.bss.Utils.sumArray(probabilities);
		probabilities[K - 1] += fallThrough;

		return probabilities;
	}// END: stickBreak

	private double rBeta() throws MathException {
		return beta.inverseCumulativeProbability(MathUtils.nextDouble());
	}// END: rbeta

	public ArrayList<Integer> chinRest(int N, double gamma) {

		ArrayList<Integer> mapping = new ArrayList<Integer>();

		mapping.add(0, 0);
		int nextFree = 1;

		for (int i = 1; i < N; i++) {

			double u1 = Math.random();
			if (u1 < (gamma / (gamma + i))) {

				// sit at new table
				mapping.add(i, nextFree);
				nextFree++;

			} else {

				// choose existing table with weights by number of customers
				int samplePos = -Integer.MAX_VALUE;
				double cumProb = 0.0;
				double u2 = Math.random();
				double size = (double) mapping.size();

				for (int j : mapping) {

					cumProb += 1 / size;

					if (u2 < cumProb) {
						samplePos = j;
						break;
					}

				}

				mapping.add(i, samplePos);

			}// END: u1 check

		}// END: i loop

		return mapping;
	}// END: chinProc

	@Override
	public double getCoercableParameter() {
		// TODO Auto-generated method stub
		return 0;
	}

	@Override
	public void setCoercableParameter(double value) {
		// TODO Auto-generated method stub

	}

	@Override
	public double getRawParameter() {
		// TODO Auto-generated method stub
		return 0;
	}

	@Override
	public String getPerformanceSuggestion() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public String getOperatorName() {
		return zParameter.getParameterName();
	}

	public static void main(String args[]) {

	}// END: main

}// END: class
