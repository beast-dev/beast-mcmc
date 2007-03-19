package dr.inference.model;

import dr.inference.parallel.MPIServices;
import dr.xml.*;

import java.util.ArrayList;
import java.util.List;

/**
 * Created by IntelliJ IDEA.
 * User: msuchard
 * Date: Jan 22, 2007
 * Time: 1:56:07 PM
 * To change this template use File | Settings | File Templates.
 */
public class ParallelCompoundLikelihood extends CompoundLikelihood {

	public static final String PARALLEL_COMPOUND_LIKELIHOOD = "parallelCompoundLikelihood";
	public static final String LOCAL_CHECK = "doLocalCheck";
	public static final String RUN_PARALLEL = "doInParallel";

	public ParallelCompoundLikelihood(boolean doParallel, boolean checkLocal) {
		super();
		this.doParallel = doParallel;
		this.checkLocal = checkLocal;
	}

	public void addLikelihood(Likelihood likelihood) {
		super.addLikelihood(likelihood);
		// todo make sure that link to a remote node exists
	}


	private boolean doParallel = true;
	private boolean checkLocal = false;

	public double getLogLikelihood() {

		double logLikelihood = 0;
		if (doParallel) {

			logLikelihood = getLogLikelihoodRemote();
			if (checkLocal) {
				super.makeDirty();
				double logLikelihoodLocal = super.getLogLikelihood();
				System.err.printf("Local: %5.4f  Remote: %5.4f\n", logLikelihoodLocal, logLikelihood);
			}

		} else
			logLikelihood = super.getLogLikelihood();

		return logLikelihood;
	}

	private double getLogLikelihoodRemote() {
		double logLikelihood = 0.0;

		final int N = getLikelihoodCount();

		List<ParallelLikelihood> likelihoodsDistributed = new ArrayList<ParallelLikelihood>();
		List<Integer> processorList = new ArrayList<Integer>();

		for (int i = 0; i < N; i++) {
			ParallelLikelihood likelihood = (ParallelLikelihood) getLikelihood(i);
			if (!likelihood.getLikelihoodKnown()) {
				//    if (true) {
				final int processor = i + 1;
//				MPIServices.requestLikelihood(processor);
//				((AbstractModel) getLikelihood(i).getModel()).sendState(processor);
				likelihoodsDistributed.add(likelihood);
				processorList.add(processor);
			} else {
				logLikelihood += likelihood.getLogLikelihood();
			}
		}

		final int size = likelihoodsDistributed.size();

		if (size == 1) { // only one, so do locally

			logLikelihood += likelihoodsDistributed.get(0).getLogLikelihood();

		} else if (size > 1) {

			// Distribute calculations
			int index = 0;
			for (ParallelLikelihood likelihood : likelihoodsDistributed) {
				int processor = processorList.get(index++);
				MPIServices.requestLikelihood(processor);
				((AbstractModel) likelihood.getModel()).sendState(processor);
			}

			// Implicit barrier

			// Collect calculations
			index = 0;
			for (ParallelLikelihood likelihood : likelihoodsDistributed) {
				int processor = processorList.get(index++);
				double l = MPIServices.receiveDouble(processor);
				logLikelihood += l;
				likelihood.setLikelihood(l);         // todo don't we need to set all of the submodels ????
			}

			// todo Use Gather instead of blocking loop

		}

		return logLikelihood;
	}


	public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

		public String getParserName() {
			return PARALLEL_COMPOUND_LIKELIHOOD;
		}

		//      public String[] getParserNames() { return new String[] { getParserName(), "posterior", "prior", "likelihood" }; }

		public Object parseXMLObject(XMLObject xo) throws XMLParseException {

			boolean doParallel = true;
			boolean checkLocal = false;

			if (xo.hasAttribute(LOCAL_CHECK)) {
				checkLocal = xo.getBooleanAttribute(LOCAL_CHECK);
			}

			if (xo.hasAttribute(RUN_PARALLEL)) {
				doParallel = xo.getBooleanAttribute(RUN_PARALLEL);
			}

			ParallelCompoundLikelihood compoundLikelihood = new ParallelCompoundLikelihood(doParallel, checkLocal);

			for (int i = 0; i < xo.getChildCount(); i++) {
				if (xo.getChild(i) instanceof Likelihood) {
					compoundLikelihood.addLikelihood((Likelihood) xo.getChild(i));
				} else {
					throw new XMLParseException("An element which is not a likelihood has been added to a " + PARALLEL_COMPOUND_LIKELIHOOD + " element");
				}
			}
			return compoundLikelihood;
		}

		//************************************************************************
		// AbstractXMLObjectParser implementation
		//************************************************************************

		public String getParserDescription() {
			return "A likelihood function which is simply the product of its component likelihood functions.";
		}

		public XMLSyntaxRule[] getSyntaxRules() {
			return rules;
		}

		private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
				new ElementRule(Likelihood.class, 1, Integer.MAX_VALUE),
		};

		public Class getReturnType() {
			return ParallelCompoundLikelihood.class;
		}
	};
}
