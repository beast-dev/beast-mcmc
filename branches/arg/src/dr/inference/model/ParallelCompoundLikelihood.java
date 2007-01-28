package dr.inference.model;

import dr.inference.parallel.MPIServices;
import dr.xml.*;

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
		if (doParallel)
			logLikelihood = getLogLikelihoodRemote();
		else
			logLikelihood = super.getLogLikelihood();

		if (checkLocal && doParallel) {
			double logLikelihoodLocal = super.getLogLikelihood();
			System.err.printf("Local: %5.4f  Remote: %5.4f\n",logLikelihoodLocal,logLikelihood);
		}

		return logLikelihood;
	}



    private double getLogLikelihoodRemote() {
        double logLikelihood = 0.0;

        final int N = getLikelihoodCount();

        for (int i = 0; i < N; i++) {

           MPIServices.requestLikelihood(i+1);
	       ((AbstractModel)getLikelihood(i).getModel()).sendState(i+1);
        }

	    // Implicit barrier

	    for (int i=0; i<N; i++) {
	        double l = MPIServices.receiveDouble(i+1);
	        if (l == Double.NEGATIVE_INFINITY) return Double.NEGATIVE_INFINITY;

            logLikelihood += l;
        }

	    // todo Use Reduce instead of blocking loop



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

            ParallelCompoundLikelihood compoundLikelihood = new ParallelCompoundLikelihood(doParallel,checkLocal);

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
