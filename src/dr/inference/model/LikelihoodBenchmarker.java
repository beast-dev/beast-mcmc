package dr.inference.model;

import dr.inference.model.Likelihood;
import dr.xml.*;

import java.util.List;
import java.util.ArrayList;
import java.util.logging.Logger;

/**
 * @author Andrew Rambaut
 * @version $Id$
 */
public class LikelihoodBenchmarker {

    public LikelihoodBenchmarker(List<Likelihood> likelihoods, int iterationCount) {
        for (Likelihood likelihood : likelihoods) {
            long startTime = System.nanoTime();

            for (int i = 0; i < iterationCount; i++) {
                likelihood.makeDirty();
                likelihood.getLogLikelihood();
            }

            long endTime = System.nanoTime();

            Logger.getLogger("dr.evomodel.beagle").info(
                    "Benchmark " + likelihood.getId() + "(" + likelihood.getClass().getName() + "): " +
                            (endTime - startTime) + "ns");

        }
    }

    public static AbstractXMLObjectParser PARSER =  new AbstractXMLObjectParser() {

        public static final String BENCHMARKER = "benchmarker";

        public String getParserName() {
            return BENCHMARKER;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {
            int iterationCount = 1000;

            if (xo.hasAttribute("iterationCount")) {
                iterationCount = xo.getIntegerAttribute("iterationCount");
            }

            List<Likelihood> likelihoods = new ArrayList<Likelihood>();

            for (int i = 0; i < xo.getChildCount(); i++) {
                Object xco = xo.getChild(i);
                if (xco instanceof Likelihood) {
                    likelihoods.add((Likelihood)xco);
                }
            }

            if (likelihoods.size() == 0) {
                throw new XMLParseException("No likelihoods for benchmarking");
            }

            return new LikelihoodBenchmarker(likelihoods, iterationCount);
        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public String getParserDescription() {
            return "This element runs a benchmark on a series of likelihood calculators.";
        }

        public Class getReturnType() {
            return Likelihood.class;
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private final XMLSyntaxRule[] rules = {
                AttributeRule.newIntegerRule("iterationCount"),
                new ElementRule(Likelihood.class, 1, Integer.MAX_VALUE)
        };
    };
}
