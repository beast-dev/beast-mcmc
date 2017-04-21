package dr.inferencexml.model;

import dr.inference.model.LFMLoadingsPotentialDerivative;
import dr.inference.model.LatentFactorModel;
import dr.xml.*;

    /**
     * @author Max Tolkoff
     */
    public class LFMLoadingsPotentialDerivativeParser extends AbstractXMLObjectParser {
        public final static String LFM_LOADINGS_DERIVATIVE = "LFMLoadingsPotentialDerivative";

        @Override
        public String getParserName() {
            return LFM_LOADINGS_DERIVATIVE;
        }

        @Override
        public Object parseXMLObject(XMLObject xo) throws XMLParseException {
            LatentFactorModel lfm = (LatentFactorModel) xo.getChild(LatentFactorModel.class);

            return new LFMLoadingsPotentialDerivative(lfm);
        }

        @Override
        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private final XMLSyntaxRule[] rules = {
                new ElementRule(LatentFactorModel.class),
        };

        @Override
        public String getParserDescription() {
            return null;
        }

        @Override
        public Class getReturnType() {
            return LFMLoadingsPotentialDerivative.class;
        }
    }

