package dr.evoxml;

import dr.xml.*;
import dr.inference.model.Parameter;
import dr.evomodel.continuous.MultivariateDiffusionModel;
import dr.evomodel.continuous.BifractionalDiffusionModel;

/**
 * @author Marc Suchard
 */

public class BifractionalDiffusionModelParser extends AbstractXMLObjectParser {

    public static final String BIFRACTIONAL_DIFFUSION_PROCESS = "bifractionalDiffusionModel";
    public static final String ALPHA_PARAMETER = "alpha";
    public static final String BETA_PARAMETER = "beta";

    public String getParserName() {
            return BIFRACTIONAL_DIFFUSION_PROCESS;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            XMLObject cxo = xo.getChild(ALPHA_PARAMETER);
            Parameter alpha = (Parameter) cxo.getChild(Parameter.class);

            cxo = xo.getChild(BETA_PARAMETER);
            Parameter beta = (Parameter) cxo.getChild(Parameter.class);

            return new BifractionalDiffusionModel(alpha, beta);
        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public String getParserDescription() {
            return "Describes a bivariate diffusion process using a bifractional random walk";
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private final XMLSyntaxRule[] rules = {
                new ElementRule(ALPHA_PARAMETER,
                        new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
                new ElementRule(BETA_PARAMETER,
                        new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
        };

        public Class getReturnType() {
            return MultivariateDiffusionModel.class;
        }
}
