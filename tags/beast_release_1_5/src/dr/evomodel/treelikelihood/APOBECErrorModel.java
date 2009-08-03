package dr.evomodel.treelikelihood;

import dr.inference.model.Parameter;
import dr.inference.model.Statistic;
import dr.xml.*;

import java.util.logging.Logger;

/**
 * @author Andrew Rambaut
 * @version $Id$
 */
public class APOBECErrorModel extends TipPartialsModel {
    public enum APOBECType {
        ALL("all"),
        BOTH("both"),
        HA3G("hA3G"),
        HA3F("hA3F");


        APOBECType(String label) {
            this.label = label;
        }

        public String toString() {
            return label;
        }

        final String label;
    }

    public static final String APOBEC_ERROR_MODEL = "APOBECErrorModel";
    public static final String HYPERMUTATION_RATE = "hypermutationRate";
    public static final String HYPERMUTATION_INDICATORS = "hypermutationIndicators";

    public APOBECErrorModel(APOBECType type, Parameter hypermutationRateParameter, Parameter hypermuationIndicatorParameter) {
        super(APOBEC_ERROR_MODEL, null, null);

        this.type = type;

        this.hypermutationRateParameter = hypermutationRateParameter;
        addVariable(this.hypermutationRateParameter);


        this.hypermuationIndicatorParameter = hypermuationIndicatorParameter;
        addVariable(this.hypermuationIndicatorParameter);

        addStatistic(new TaxonHypermutatedStatistic());
    }

    public void getTipPartials(int nodeIndex, double[] partials) {
        int[] states = this.states[nodeIndex];

        if (hypermuationIndicatorParameter.getParameterValue(nodeIndex) > 0.0) {
            double rate = hypermutationRateParameter.getParameterValue(0);

            int k = 0;
            int nextState;
            for (int j = 0; j < patternCount; j++) {

                switch (states[j]) {
                    case 0: // is an A
                        double pMutated = 0.0;
                        if (j < patternCount - 1) {
                            nextState = states[j+1];

                            if (    (type == APOBECType.ALL) ||
                                    (type == APOBECType.HA3G && nextState == 2) || // is a G
                                    (type == APOBECType.HA3F && nextState == 0) || // is an A
                                    (type == APOBECType.BOTH && (nextState == 2 || nextState == 0))
                                    ) {
                                pMutated = rate;
                            }
                        }
                        partials[k] = 1.0 - pMutated;
                        partials[k + 1] = 0.0;
                        partials[k + 2] = pMutated;
                        partials[k + 3] = 0.0;

                        break;
                    case 1: // is an C
                        partials[k] = 0.0;
                        partials[k + 1] = 1.0;
                        partials[k + 2] = 0.0;
                        partials[k + 3] = 0.0;
                        break;
                    case 2: // is an G
                        partials[k] = 0.0;
                        partials[k + 1] = 0.0;
                        partials[k + 2] = 1.0;
                        partials[k + 3] = 0.0;
                        break;
                    case 3: // is an T
                        partials[k] = 0.0;
                        partials[k + 1] = 0.0;
                        partials[k + 2] = 0.0;
                        partials[k + 3] = 1.0;
                        break;
                    default: // is an ambiguity
                        partials[k] = 1.0;
                        partials[k + 1] = 1.0;
                        partials[k + 2] = 1.0;
                        partials[k + 3] = 1.0;
                }

                k += stateCount;
            }
        } else {

            int k = 0;
            for (int j = 0; j < patternCount; j++) {
                switch (states[j]) {
                    case 0: // is an A
                        partials[k] = 1.0;
                        partials[k + 1] = 0.0;
                        partials[k + 2] = 0.0;
                        partials[k + 3] = 0.0;
                        break;
                    case 1: // is an C
                        partials[k] = 0.0;
                        partials[k + 1] = 1.0;
                        partials[k + 2] = 0.0;
                        partials[k + 3] = 0.0;
                        break;
                    case 2: // is an G
                        partials[k] = 0.0;
                        partials[k + 1] = 0.0;
                        partials[k + 2] = 1.0;
                        partials[k + 3] = 0.0;
                        break;
                    case 3: // is an T
                        partials[k] = 0.0;
                        partials[k + 1] = 0.0;
                        partials[k + 2] = 0.0;
                        partials[k + 3] = 1.0;
                        break;
                    default: // is an ambiguity
                        partials[k] = 1.0;
                        partials[k + 1] = 1.0;
                        partials[k + 2] = 1.0;
                        partials[k + 3] = 1.0;
                }
                k += stateCount;
            }
        }

    }

    public class TaxonHypermutatedStatistic extends Statistic.Abstract {

        public TaxonHypermutatedStatistic() {
            super("hypermutated");
        }

        public int getDimension() {
            return hypermuationIndicatorParameter.getDimension();
        }

        public String getDimensionName(int dim) {
            return taxonMap.get(dim);
        }

        public double getStatisticValue(int i) {
            return hypermuationIndicatorParameter.getParameterValue(i);
        }

    }

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public String getParserName() { return APOBEC_ERROR_MODEL; }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            APOBECType type = APOBECType.HA3G;

            if (xo.hasAttribute("type")) {
                if (xo.getStringAttribute("type").equalsIgnoreCase("all")) {
                    type = APOBECType.ALL;
                } else if (xo.getStringAttribute("type").equalsIgnoreCase("both")) {
                    type = APOBECType.BOTH;
                } else if (xo.getStringAttribute("type").equalsIgnoreCase("hA3F")) {
                    type = APOBECType.HA3F;
                } else if (!xo.getStringAttribute("type").equalsIgnoreCase("hA3G")) {
                    throw new XMLParseException("unrecognized option for attribute, 'type': " + xo.getStringAttribute("type"));
                }
            }

            Parameter hypermutationRateParameter = null;
            if (xo.hasChildNamed(HYPERMUTATION_RATE)) {
                hypermutationRateParameter = (Parameter)xo.getElementFirstChild(HYPERMUTATION_RATE);
            }

            Parameter hypermuationIndicatorParameter = null;
            if (xo.hasChildNamed(HYPERMUTATION_INDICATORS)) {
                hypermuationIndicatorParameter = (Parameter)xo.getElementFirstChild(HYPERMUTATION_INDICATORS);
            }

            APOBECErrorModel errorModel =  new APOBECErrorModel(
                    type, hypermutationRateParameter, hypermuationIndicatorParameter);

            Logger.getLogger("dr.evomodel").info("Using APOBEC error model, assuming APOBEC " + type.name());

            return errorModel;
        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public String getParserDescription() {
            return
                    "This element returns a model that allows for APOBEC-type RNA editing.";
        }

        public Class getReturnType() { return APOBECErrorModel.class; }

        public XMLSyntaxRule[] getSyntaxRules() { return rules; }

        private XMLSyntaxRule[] rules = new XMLSyntaxRule[] {
                AttributeRule.newStringRule("type", true),
                new ElementRule(HYPERMUTATION_RATE, Parameter.class, "The hypermutation rate per target site per sequence"),
                new ElementRule(HYPERMUTATION_INDICATORS, Parameter.class, "A binary indicator of whether the sequence is hypermutated"),
        };
    };

    private final APOBECType type;
    private final Parameter hypermutationRateParameter;
    private final Parameter hypermuationIndicatorParameter;
}