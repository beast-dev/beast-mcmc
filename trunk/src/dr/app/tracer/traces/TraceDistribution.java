package dr.app.tracer.traces;

import dr.stats.DiscreteStatistics;
import dr.util.HeapSort;

/**
 * A class that stores the distribution statistics for a trace
 *
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @version $Id: TraceDistribution.java,v 1.1.1.2 2006/04/25 23:00:09 rambaut Exp $
 */
public class TraceDistribution {

	public TraceDistribution(double[] values, int stepSize) {
		analyseDistribution(values, stepSize);
	}

	public TraceDistribution(double[] values, int stepSize, double ESS) {
		analyseDistribution(values, stepSize);
		this.ESS = ESS;
	}

	public boolean isValid() { return isValid; }

	public double getMean() { return mean; }
	public double getMedian() { return median; }
	public double getLowerHPD() { return hpdLower; }
	public double getUpperHPD() { return hpdUpper; }
	public double getLowerCPD() { return cpdLower; }
	public double getUpperCPD() { return cpdUpper; }
	public double getESS() { return ESS; }

	public double getMinimum() { return minimum; }
	public double getMaximum() { return maximum; }
	/**
	 * Analyze trace
	 */
	private void analyseDistribution(double[] values, int stepSize) {

		mean = DiscreteStatistics.mean(values);

		minimum = Double.POSITIVE_INFINITY;
		maximum = Double.NEGATIVE_INFINITY;

		for (int i = 0; i < values.length; i++) {
			if (values[i] < minimum) minimum = values[i];
			if (values[i] > maximum) maximum = values[i];
		}

		if (maximum == minimum) {
			isValid = false;
			return;
		}

		int[] indices = new int[values.length];
		HeapSort.sort(values, indices);
		median = DiscreteStatistics.quantile(0.5, values, indices);
		cpdLower = DiscreteStatistics.quantile(0.025, values, indices);
		cpdUpper = DiscreteStatistics.quantile(0.975, values, indices);
		calculateHPDInterval(0.95, values, indices);
		ESS = values.length;

		isValid = true;
	}

	/**
	 * @param proportion the proportion of probability mass oncluded within interval.
	 */
	private void calculateHPDInterval(double proportion, double[] array, int[] indices) {

		double minRange = Double.MAX_VALUE;
		int hpdIndex = 0;

		int diff = (int)Math.round(proportion * (double)array.length);
		for (int i =0; i <= (array.length - diff); i++) {
			double minValue = array[indices[i]];
			double maxValue = array[indices[i+diff-1]];
			double range = Math.abs(maxValue - minValue);
			if (range < minRange) {
				minRange = range;
				hpdIndex = i;
			}
		}
		hpdLower = array[indices[hpdIndex]];
		hpdUpper = array[indices[hpdIndex+diff-1]];
	}

	//************************************************************************
	// private methods
	//************************************************************************

	protected boolean isValid = false;

	protected double minimum, maximum;
	protected double mean, median;
	protected double cpdLower, cpdUpper, hpdLower, hpdUpper;
	protected double ESS;
}