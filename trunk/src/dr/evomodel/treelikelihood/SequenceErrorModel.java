package dr.evomodel.treelikelihood;

import dr.evolution.util.TaxonList;
import dr.inference.model.Parameter;
import dr.xml.*;

import java.util.logging.Logger;

/**
 * @author Andrew Rambaut
 * @version $Id$
 */
public class SequenceErrorModel extends TipPartialsModel {
    public enum ErrorType {
        TYPE_1_TRANSITIONS("type1Transitions"),
        TYPE_2_TRANSITIONS("type2Transitions"),
        TRANSITIONS_ONLY("transitionsOnly"),
        ALL_SUBSTITUTIONS("allSubstitutions");


        ErrorType(String label) {
            this.label = label;
        }

        public String toString() {
            return label;
        }

        final String label;
    }

    public static final String SEQUENCE_ERROR_MODEL = "sequenceErrorModel";
    public static final String BASE_ERROR_RATE = "baseErrorRate";
    public static final String AGE_RELATED_RATE = "ageRelatedErrorRate";

    public static final String EXCLUDE = "exclude";
    public static final String INCLUDE = "include";

    public static final String TYPE = "type";

    public SequenceErrorModel(TaxonList includeTaxa, TaxonList excludeTaxa,
                              ErrorType errorType, Parameter baseErrorRateParameter, Parameter ageRelatedErrorRateParameter) {
        super(SEQUENCE_ERROR_MODEL, includeTaxa, excludeTaxa);

        this.errorType = errorType;

        if (baseErrorRateParameter != null) {
            this.baseErrorRateParameter = baseErrorRateParameter;
            addParameter(this.baseErrorRateParameter);
        } else {
            this.baseErrorRateParameter = null;
        }

        if (ageRelatedErrorRateParameter != null) {
            this.ageRelatedErrorRateParameter = ageRelatedErrorRateParameter;
            addParameter(ageRelatedErrorRateParameter);
        } else {
            this.ageRelatedErrorRateParameter = null;
        }
    }

    public void getTipPartials(int nodeIndex, double[] partials) {
        int[] states = this.states[nodeIndex];

        double pUndamaged = 1.0;
        double pDamagedTS = 0.0;
        double pDamagedTV = 0.0;

        if (!excluded[nodeIndex]) {
            if (baseErrorRateParameter != null) {
                pUndamaged = pUndamaged - baseErrorRateParameter.getParameterValue(0);
            }

            if (ageRelatedErrorRateParameter != null) {
                double rate = ageRelatedErrorRateParameter.getParameterValue(0);
                double age = tree.getNodeHeight(tree.getExternalNode(nodeIndex));
                pUndamaged *= Math.exp(-rate * age);
            }


            if (errorType == ErrorType.ALL_SUBSTITUTIONS) {
                pDamagedTS = (1.0 - pUndamaged) / 3.0;
                pDamagedTV = pDamagedTS;

            } else if (errorType == ErrorType.TRANSITIONS_ONLY) {
                pDamagedTS = 1.0 - pUndamaged;
                pDamagedTV = 0.0;
            } else {
                throw new IllegalArgumentException("only TRANSITIONS_ONLY and ALL_SUBSTITUTIONS are supported");
            }

        }
        
        int k = 0;
        for (int j = 0; j < patternCount; j++) {
            switch (states[j]) {
                case 0: // is an A
                    partials[k] = pUndamaged;
                    partials[k + 1] = pDamagedTV;
                    partials[k + 2] = pDamagedTS;
                    partials[k + 3] = pDamagedTV;
                    break;
                case 1: // is an C
                    partials[k] = pDamagedTV;
                    partials[k + 1] = pUndamaged;
                    partials[k + 2] = pDamagedTV;
                    partials[k + 3] = pDamagedTS;
                    break;
                case 2: // is an G
                    partials[k] = pDamagedTS;
                    partials[k + 1] = pDamagedTV;
                    partials[k + 2] = pUndamaged;
                    partials[k + 3] = pDamagedTV;
                    break;
                case 3: // is an T
                    partials[k] = pDamagedTV;
                    partials[k + 1] = pDamagedTS;
                    partials[k + 2] = pDamagedTV;
                    partials[k + 3] = pUndamaged;
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

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public String getParserName() { return SEQUENCE_ERROR_MODEL; }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            ErrorType errorType = ErrorType.ALL_SUBSTITUTIONS;

            if (xo.hasAttribute(TYPE)) {
                if (xo.getStringAttribute(TYPE).equalsIgnoreCase("transitions")) {
                    errorType = ErrorType.TRANSITIONS_ONLY;
                } else if (!xo.getStringAttribute(TYPE).equalsIgnoreCase("all")) {
                    throw new XMLParseException("unrecognized option for attribute, 'type': " + xo.getStringAttribute(TYPE));
                }
            }

            Parameter baseDamageRateParameter = null;
            if (xo.hasChildNamed(BASE_ERROR_RATE)) {
                baseDamageRateParameter = (Parameter)xo.getElementFirstChild(BASE_ERROR_RATE);
            }

            Parameter ageRelatedRateParameter = null;
            if (xo.hasChildNamed(AGE_RELATED_RATE)) {
                ageRelatedRateParameter = (Parameter)xo.getElementFirstChild(AGE_RELATED_RATE);
            }

            if (baseDamageRateParameter == null && ageRelatedRateParameter == null) {
                throw new XMLParseException("You must specify one or other or both of " +
                        BASE_ERROR_RATE + " and " + AGE_RELATED_RATE + " parameters");
            }

            TaxonList includeTaxa = null;
            TaxonList excludeTaxa = null;

            if (xo.hasChildNamed(INCLUDE)) {
                includeTaxa = (TaxonList)xo.getElementFirstChild(INCLUDE);
            }

            if (xo.hasChildNamed(EXCLUDE)) {
                excludeTaxa = (TaxonList)xo.getElementFirstChild(EXCLUDE);
            }

            SequenceErrorModel aDNADamageModel =  new SequenceErrorModel(includeTaxa, excludeTaxa,
                    errorType, baseDamageRateParameter, ageRelatedRateParameter);

            Logger.getLogger("dr.evomodel").info("Using sequence error model, assuming errors cause " +
                    (errorType == ErrorType.TRANSITIONS_ONLY ? "transitions only." : "any substitution."));

            return aDNADamageModel;
        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public String getParserDescription() {
            return "This element returns a model that allows for post-mortem DNA damage.";
        }

        public Class getReturnType() { return SequenceErrorModel.class; }

        public XMLSyntaxRule[] getSyntaxRules() { return rules; }

        private final XMLSyntaxRule[] rules = {
                AttributeRule.newStringRule(TYPE, true),
                new ElementRule(BASE_ERROR_RATE, Parameter.class, "The base error rate per site per sequence", true),
                new ElementRule(AGE_RELATED_RATE, Parameter.class, "The error rate per site per unit time", true),
                new XORRule(
                        new ElementRule(INCLUDE, TaxonList.class, "A set of taxa to which to apply the damage model to"),
                        new ElementRule(EXCLUDE, TaxonList.class, "A set of taxa to which to not apply the damage model to")
                        , true)
        };
    };

    private final ErrorType errorType;
    private final Parameter baseErrorRateParameter;
    private final Parameter ageRelatedErrorRateParameter;
}