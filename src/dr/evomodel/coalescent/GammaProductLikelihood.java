package dr.evomodel.coalescent;

import dr.evomodel.tree.TreeModel;
import dr.inference.model.Likelihood;
import dr.inference.model.Parameter;
import dr.math.distributions.GammaDistribution;

/**
 * Calculates a product of gamma distributions.
 *
 * @author Guy Baele
 */

public class GammaProductLikelihood extends Likelihood.Abstract {
	
	public final static boolean REDUCE_TO_EXPONENTIAL = true;
	
	private TreeModel treeModel;
	private double popSize;
	private double[] means;
	private double[] variances;
	
	public GammaProductLikelihood(TreeModel treeModel, double popSize, double[] means, double[] variances) {
		super(treeModel);
		this.treeModel = treeModel;
		this.popSize = popSize;
		this.means = means;
		this.variances = variances;
	}

	public double calculateLogLikelihood() {
		
		double logPDF = 0.0;
		
		//means and variances are probably in the reverse order
		System.err.println("\nProvided empirical means and variances: ");
		for (int i = 0; i < means.length; i++) {
			System.err.println(means[i] + "   " + variances[i]);
		}
		CoalescentTreeIntervalStatistic ctis = new CoalescentTreeIntervalStatistic(treeModel);
		System.err.println("\nDimension = " + ctis.getDimension() + "\nLineage info: ");
		/*for (int i = 0; i < ctis.getDimension(); i++) {
			System.err.println(ctis.getLineageCount(i));
		}*/
		for (int i = ctis.getDimension()-1; i >= 0; i--) {
			System.err.println(ctis.getLineageCount(i));
		}
		System.err.println("\nStatistic values: ");
		for (int i = ctis.getDimension()-1; i >= 0; i--) {
			System.err.println(ctis.getStatisticValue(i));
		}
		
		//ignore possibility of more than 1 dimension for now
		//System.err.println("\nPopulation size = " + popSize);
		System.err.println("\nTree: " + treeModel);
		
		//calculate alpha and beta for the gamma distributions
		double[] alphas = new double[means.length];
		for (int i = 0; i < alphas.length; i++) {
			alphas[i] = means[i]*means[i]/variances[i];
		}
		double[] betas = new double[means.length];
		for (int i = 0; i < betas.length; i++) {
			betas[i] = means[i]/variances[i];
		}
		
		if (REDUCE_TO_EXPONENTIAL) {
			for (int i = 0; i < alphas.length; i++) {
				alphas[i] = 1.0;
			}
			for (int i = 0; i < betas.length; i++) {
				betas[i] = 1.0/popSize;
			}
		}
		
		int indicator = 0;
		for (int i = ctis.getDimension()-1; i > 0; i--) {
			System.err.println("\nInterval " + i);
			if (i == ctis.getDimension()-1) {
				//coalescent event: gamma density
				System.err.println("Coalescent event at root");
				System.err.println("Interval length = " + ctis.getStatisticValue(i));
				System.err.println("Lineage count = " + ctis.getLineageCount(i));
				int combinations = (int)ctis.getLineageCount(i)*((int)ctis.getLineageCount(i)-1)/2;
				System.err.println("Combinations = " + combinations);
				logPDF += GammaDistribution.logPdf(ctis.getStatisticValue(i), alphas[indicator], 1.0/(combinations*betas[indicator]));
				System.err.println(GammaDistribution.logPdf(ctis.getStatisticValue(i), alphas[indicator], 1.0/(combinations*betas[indicator])));
				if (indicator < (means.length-1)) {
					indicator++;
				}
			} else if (ctis.getLineageCount(i) - ctis.getLineageCount(i+1) > 0) {
				//coalescent event: gamma density
				System.err.println("Coalescent event");
				System.err.println("Interval length = " + ctis.getStatisticValue(i));
				System.err.println("Lineage count = " + ctis.getLineageCount(i));
				int combinations = (int)ctis.getLineageCount(i)*((int)ctis.getLineageCount(i)-1)/2;
				System.err.println("Combinations = " + combinations);
				logPDF += GammaDistribution.logPdf(ctis.getStatisticValue(i), alphas[indicator], 1.0/(combinations*betas[indicator]));
				System.err.println(GammaDistribution.logPdf(ctis.getStatisticValue(i), alphas[indicator], 1.0/(combinations*betas[indicator])));
				if (indicator < (means.length-1)) {
					indicator++;
				}
			} else {
				//new sample: gamma tail probability
				System.err.println("New sample");
				System.err.println("Interval length = " + ctis.getStatisticValue(i));
				System.err.println("Lineage count = " + ctis.getLineageCount(i));
				int combinations = (int)ctis.getLineageCount(i)*((int)ctis.getLineageCount(i)-1)/2;
				System.err.println("Combinations = " + combinations);
				logPDF += Math.log(1-GammaDistribution.cdf(ctis.getStatisticValue(i), alphas[indicator], 1.0/(combinations*betas[indicator])));
				System.err.println(Math.log(1-GammaDistribution.cdf(ctis.getStatisticValue(i), alphas[indicator], 1.0/(combinations*betas[indicator]))));
			}
		}
		
		System.err.println("\nlogPDF = " + logPDF);
		System.exit(0);
		
		return logPDF;
	}

}
