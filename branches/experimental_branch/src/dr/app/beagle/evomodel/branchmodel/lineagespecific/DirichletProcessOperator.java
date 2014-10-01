package dr.app.beagle.evomodel.branchmodel.lineagespecific;

import org.apache.commons.math.MathException;

import dr.inference.model.Parameter;
import dr.inference.operators.AbstractCoercableOperator;
import dr.inference.operators.CoercionMode;
import dr.inference.operators.OperatorFailedException;
import dr.inference.operators.SimpleMCMCOperator;
import dr.math.MathUtils;

public class DirichletProcessOperator extends SimpleMCMCOperator
//AbstractCoercableOperator
{

//	http://xiaodong-yu.blogspot.be/2009/09/gibbs-sampling-for-dp-mixtures.html
//	http://www.arbylon.net/resources.html

	private DirichletProcessPrior dpp;
	
	private Parameter zParameter;
	private int realizationCount;
	private int uniqueRealizationCount;
//	private double intensity;

	public DirichletProcessOperator(DirichletProcessPrior dpp, Parameter zParameter, int uniqueRealizationCount, double weight) {
		this(dpp, zParameter, null,
				uniqueRealizationCount, weight);
	}// END: Constructor

	public DirichletProcessOperator(DirichletProcessPrior dpp, Parameter zParameter, CoercionMode mode, 
			int uniqueRealizationCount, double weight) {

//		super(mode);

		this.dpp = dpp;
		
		this.zParameter = zParameter;
//		this.intensity = intensity;
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
	
	
//	private void doOperate() throws MathException {
//		
//		int index = MathUtils.nextInt(realizationCount);
//		
//		// compute configuration probabilities
//		double[] probs = new double[uniqueRealizationCount];
//
//		for (int i = 0; i < realizationCount; i++) {
//			probs[(int) zParameter.getParameterValue(i)]++;
//		}// END: i loop
//		
//		for (int i = 0; i < uniqueRealizationCount; i++) {
//			probs[i] = probs[i] / realizationCount;
//		}// END: i loop
//		
////		dr.app.bss.Utils.printArray(probs);
//		
//		
//		int categoryIndex = MathUtils.randomChoicePDF(probs);
//		
//		System.out.println("z[" + index + "]=" + categoryIndex);
//		
//		zParameter.setParameterValue(index, categoryIndex);
//		
//		
//	}//END: doOperate

	private void doOperate() throws MathException {

		int index = MathUtils.nextInt(realizationCount);
		
//		for (int index = 0; index < realizationCount; index++) {

       double intensity = dpp.getGamma();     
		
//       System.out.println("intensity: " + intensity);
       
			int unoccupied = uniqueRealizationCount;
			int[] occupancy = new int[uniqueRealizationCount];

			for (int i = 0; i < realizationCount; i++) {

				if (i != index) {

					int j = (int) zParameter.getParameterValue(i);
					occupancy[j]++;

					if (occupancy[j] == 1) {

						unoccupied--;

					}// END: first check

				}// END: i check

			}// END: i loop

//			dr.app.bss.Utils.printArray(occupancy);
//			System.out.println(unoccupied);
			
			double p1 = intensity / ((realizationCount - 1 + intensity) * unoccupied);
			double p = 0;
			double[] probs = new double[uniqueRealizationCount];

			for (int i = 0; i < uniqueRealizationCount; i++) {

				if (occupancy[i] == 0) {

					p = p1;

				} else {

					p = occupancy[i] / (realizationCount - 1 + intensity);

				}

				probs[i] = Math.log(p);

			}// END: i loop
			
            for (int i = 0; i < uniqueRealizationCount; i++) {
                zParameter.setParameterValue(index, i);
//                probs[i] +=  dpp.getLogLikelihood();
            }
			
			this.rescale(probs);
			this.exp(probs);

			int categoryIndex = MathUtils.randomChoicePDF(probs);
			
			zParameter.setParameterValue(index, categoryIndex);

//			 System.out.println("z[" + index + "]=" + categoryIndex);
			// System.out.println(zParameter.getParameterValue(index));

			// System.exit(0);

//		}// END: index loop

	}// END: doOperate
	
    private void exp(double[] logX) {
        for (int i = 0; i < logX.length; ++i) {
            logX[i] = Math.exp(logX[i]);
        }
    }

    private void rescale(double[] logX) {
        double max = this.max(logX);
        for (int i = 0; i < logX.length; ++i) {
            logX[i] -= max;
        }
    }

    private double max(double[] x) {
        double max = x[0];
        for (double xi : x) {
            if (xi > max) {
                max = xi;
            }
        }
        return max;
    }
	
//	@Override
//	public double getCoercableParameter() {
//		// TODO Auto-generated method stub
//		return 0;
//	}
//
//	@Override
//	public void setCoercableParameter(double value) {
//		// TODO Auto-generated method stub
//
//	}
//
//	@Override
//	public double getRawParameter() {
//		// TODO Auto-generated method stub
//		return 0;
//	}

	@Override
	public String getPerformanceSuggestion() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public String getOperatorName() {
		return zParameter.getParameterName();
	}

}// END: class
