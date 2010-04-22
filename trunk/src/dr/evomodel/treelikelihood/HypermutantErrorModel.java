package dr.evomodel.treelikelihood;

import dr.evolution.datatype.Nucleotides;
import dr.inference.model.Parameter;
import dr.inference.model.Statistic;
import dr.xml.*;

import java.util.logging.Logger;

/**
 * @author Andrew Rambaut
 * @version $Id$
 */
public class HypermutantErrorModel extends TipPartialsModel {

    public static final String HYPERMUTANT_ERROR_MODEL = "hypermutantErrorModel";
    public static final String HYPERMUTATION_RATE = "hypermutationRate";
    public static final String HYPERMUTATION_INDICATORS = "hypermutationIndicators";

    public HypermutantErrorModel(Parameter hypermutationRateParameter, Parameter hypermuationIndicatorParameter) {
        super(HYPERMUTANT_ERROR_MODEL, null, null);
        this.hypermutationRateParameter = hypermutationRateParameter;
        addVariable(this.hypermutationRateParameter);

        this.hypermuationIndicatorParameter = hypermuationIndicatorParameter;

        addVariable(this.hypermuationIndicatorParameter);

        addStatistic(new TaxonHypermutatedStatistic());
    }

    protected void taxaChanged() {
        if (hypermuationIndicatorParameter.getDimension() <= 1) {
            this.hypermuationIndicatorParameter.setDimension(tree.getExternalNodeCount());
        }
    }

    public void getTipPartials(int nodeIndex, double[] partials) {
        int[] states = this.states[nodeIndex];
        boolean isHypermutated = hypermuationIndicatorParameter.getParameterValue(nodeIndex) > 0.0;

        double rate = hypermutationRateParameter.getParameterValue(0);

        int k = 0;
        for (int j = 0; j < patternCount; j++) {

            switch (states[j]) {
                case Nucleotides.A_STATE: // is an A
                    partials[k] = 1.0;
                    partials[k + 1] = 0.0;
                    partials[k + 2] = 0.0;
                    partials[k + 3] = 0.0;
                    break;
                case Nucleotides.C_STATE: // is an C
                    partials[k] = 0.0;
                    partials[k + 1] = 1.0;
                    partials[k + 2] = 0.0;
                    partials[k + 3] = 0.0;
                    break;
                case Nucleotides.G_STATE: // is an G
                    partials[k] = 0.0;
                    partials[k + 1] = 0.0;
                    partials[k + 2] = 1.0;
                    partials[k + 3] = 0.0;
                    break;
                case Nucleotides.UT_STATE: // is an T
                    partials[k] = 0.0;
                    partials[k + 1] = 0.0;
                    partials[k + 2] = 0.0;
                    partials[k + 3] = 1.0;
                    break;
                case Nucleotides.R_STATE: // is an A in a APOBEC context
                    if (isHypermutated) {
                        partials[k] = 1.0 - rate;
                        partials[k + 1] = 0.0;
                        partials[k + 2] = rate;
                        partials[k + 3] = 0.0;
                    } else {
                        partials[k] = 1.0;
                        partials[k + 1] = 0.0;
                        partials[k + 2] = 0.0;
                        partials[k + 3] = 0.0;
                    }

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

    /*
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
    */

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public String getParserName() { return HYPERMUTANT_ERROR_MODEL; }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            Parameter hypermutationRateParameter = null;
            if (xo.hasChildNamed(HYPERMUTATION_RATE)) {
                hypermutationRateParameter = (Parameter)xo.getElementFirstChild(HYPERMUTATION_RATE);
            }

            Parameter hypermuationIndicatorParameter = null;
            if (xo.hasChildNamed(HYPERMUTATION_INDICATORS)) {
                hypermuationIndicatorParameter = (Parameter)xo.getElementFirstChild(HYPERMUTATION_INDICATORS);
            }

            HypermutantErrorModel errorModel =  new HypermutantErrorModel(hypermutationRateParameter, hypermuationIndicatorParameter);

            Logger.getLogger("dr.evomodel").info("Using APOBEC error model");

            return errorModel;
        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public String getParserDescription() {
            return
                    "This element returns a model that allows for APOBEC-type RNA editing.";
        }

        public Class getReturnType() { return HypermutantErrorModel.class; }

        public XMLSyntaxRule[] getSyntaxRules() { return rules; }

        private XMLSyntaxRule[] rules = new XMLSyntaxRule[] {
                new ElementRule(HYPERMUTATION_RATE, Parameter.class, "The hypermutation rate per target site per sequence"),
                new ElementRule(HYPERMUTATION_INDICATORS, Parameter.class, "A binary indicator of whether the sequence is hypermutated"),
        };
    };

    public class TaxonHypermutatedStatistic extends Statistic.Abstract {

        public TaxonHypermutatedStatistic() {
            super("isHypermutated");
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


    private final Parameter hypermutationRateParameter;
    private final Parameter hypermuationIndicatorParameter;
}