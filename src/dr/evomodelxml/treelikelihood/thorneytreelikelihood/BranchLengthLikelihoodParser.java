package dr.evomodelxml.treelikelihood.thorneytreelikelihood;

import dr.evomodel.treelikelihood.thorneytreelikelihood.BranchLengthLikelihoodDelegate;

import dr.xml.*;

public class BranchLengthLikelihoodParser extends AbstractXMLObjectParser {

        public static final String BRANCHLENGTH_LIKELIHOOD = "branchLengthLikelihood";

        public String getParserName() {
            return BRANCHLENGTH_LIKELIHOOD;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {
            double scale = xo.getAttribute("scale",1.0);
            return new BranchLengthLikelihoodDelegate(scale);
        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public String getParserDescription() {
            return "This element provides the likelihood of observing a branchlength given a substituion rate that is provided by a branch rate model.";
        }

        public Class getReturnType() {
            return BranchLengthLikelihoodDelegate.class;
        }

        public static final XMLSyntaxRule[] rules = {
                AttributeRule.newDoubleRule("scale",true,"a scale factor to muliply by the rate such as sequence length. default is 1"),

        };

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }
    }

