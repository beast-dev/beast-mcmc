package dr.evomodel.treelikelihood;

import dr.inference.model.Likelihood;
import dr.util.Timer;
import dr.xml.*;

import java.util.logging.Logger;

/**
 * @author Marc A. Suchard
 */
public class GPUProfiler {

	public GPUProfiler(TreeLikelihood treeLikelihood) {

		this.treeLikelihood = treeLikelihood;
		this.core = (GPUGeneralLikelihoodCore) treeLikelihood.getLikelihoodCore();
		logger = Logger.getLogger("dr.evomodel.treelikelihood");
		logger.info("Constructed profiler for " + treeLikelihood.getId());
	}

	public void run(int iterations) {

		logger.info("Starting time-trial");

		Timer timer = new Timer();
		timer.start();
		for (int i = 0; i < iterations; i++) {
			treeLikelihood.updateAllNodes(); // todo Consider more realistic case when only a few nodes need updating
			treeLikelihood.calculateLogLikelihood();
		}

		logger.info("Trial time = " + timer.toString());


	}

	public void setMode(int mode) {
		core.setMode(mode);
	}

	public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

		public static final String PROFILER_NAME = "profileGPU";
		public static final String MODE = "mode";
		public static final String TIME_TRIAL = "timeTrial";
		public static final String ITERATIONS = "iterations";

		public String getParserName() {
			return PROFILER_NAME;
		}

		public Object parseXMLObject(XMLObject xo) throws XMLParseException {

			TreeLikelihood likelihood = (TreeLikelihood) xo.getChild(TreeLikelihood.class);
			GPUProfiler profiler;

			if (likelihood.getLikelihoodCore() instanceof GPUGeneralLikelihoodCore)
				profiler = new GPUProfiler(likelihood);
			else
				throw new RuntimeException("Likelihood core in '" + likelihood.getId() + "' must use a GPU");

			int mode = xo.getAttribute(MODE, 0);

			profiler.setMode(mode);

			if (xo.getAttribute(TIME_TRIAL, false)) {
				int len = xo.getAttribute(ITERATIONS, 1000);
				profiler.run(len);
			}

			return profiler;

		}

		public String getParserDescription() {
			return "This element represents a profiling unit to test GPU algorithms.";
		}

		public Class getReturnType() {
			return Likelihood.class;
		}

		public XMLSyntaxRule[] getSyntaxRules() {
			return rules;
		}

		private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{

				new ElementRule(TreeLikelihood.class),
				AttributeRule.newIntegerRule(MODE, true),
				AttributeRule.newBooleanRule(TIME_TRIAL, true),
				AttributeRule.newIntegerRule(ITERATIONS, true),

		};
	};

	private TreeLikelihood treeLikelihood;
	private GPUGeneralLikelihoodCore core;
	private Logger logger;

}
