package dr.app.beagle.tools;

import dr.app.beagle.evomodel.treelikelihood.BeagleTreeLikelihood;
import dr.inference.loggers.LogColumn;
import dr.inference.loggers.Loggable;
import dr.inference.loggers.NumberColumn;
import dr.xml.Reportable;

public class SiteLogLikelihoodLogger implements Loggable, Reportable {

	BeagleTreeLikelihood beagleTreeLikelihood;
	int patternCount;
	private SiteLogLikelihoodColumn[] columns = null;

	public SiteLogLikelihoodLogger(BeagleTreeLikelihood beagleTreeLikelihood) {
		this.beagleTreeLikelihood = beagleTreeLikelihood;

		patternCount = beagleTreeLikelihood.getPatternCount();

	}// END: Constructor

	@Override
	public LogColumn[] getColumns() {

		if (columns == null) {
			columns = new SiteLogLikelihoodColumn[patternCount];
			for (int site = 0; site < patternCount; site++) {

				columns[site] = new SiteLogLikelihoodColumn(site);

			}
		}

		return columns;
	}// END: getColumns

	private double getSiteLogLikelihood(int site) {
		double[] siteLikelihoods = beagleTreeLikelihood.getSiteLogLikelihoods();
		return siteLikelihoods[site];
	}// END: getSiteLogLikelihoods

	private class SiteLogLikelihoodColumn extends NumberColumn {

		final int site;

		public SiteLogLikelihoodColumn(int site) {
			super("SiteLogLikelihoodColumn");
			this.site = site;
		}

		@Override
		public double getDoubleValue() {
			return getSiteLogLikelihood(site);
		}

	}// END: SiteLogLikelihoodColumn class

	public String toString() {
		StringBuilder sb = new StringBuilder();
		for (int site = 0; site < patternCount; ++site) {
			if (site > 0) {
				sb.append(", ");
			}

			sb.append(getColumns()[site].getFormatted());
		}

		return sb.toString();
	}// END: toString

	@Override
	public String getReport() {
		return toString();
	}// END: getReport

}// END: class
