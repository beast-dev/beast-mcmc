package dr.inference.model;

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

    public ParallelCompoundLikelihood() {
        super();
    }

//    @Override
//    public double getLogLikelihood() {
//        return 0;
//    }

    public void addLikelihood(Likelihood likelihood) {
        super.addLikelihood(likelihood);

        // todo link with node
    }


    public double getLogLikelihood() {
        double logLikelihood = 0.0;

        // todo distribute calculations to node and then wait for all to reply.

        final int N = getLikelihoodCount();

        for (int i = 0; i < N; i++) {
            double l = getLikelihood(i).getLogLikelihood();

            // if the likelihood is zero then short cut the rest of the likelihoods
            // This means that expensive likelihoods such as TreeLikelihoods should
            // be put after cheap ones such as BooleanLikelihoods
            if (l == Double.NEGATIVE_INFINITY) return Double.NEGATIVE_INFINITY;

            logLikelihood += l;
        }

        return logLikelihood;
    }


    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return PARALLEL_COMPOUND_LIKELIHOOD;
        }

        //      public String[] getParserNames() { return new String[] { getParserName(), "posterior", "prior", "likelihood" }; }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            // todo initialize communication with parallel nodes and distribute jobs


            ParallelCompoundLikelihood compoundLikelihood = new ParallelCompoundLikelihood();

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
