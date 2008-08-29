package dr.evomodel.treelikelihood;

import dr.evolution.util.TaxonList;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.Parameter;
import dr.xml.*;

/**
 * @author Andrew Rambaut
 * @version $Id$
 */
public class SequenceErrorModel extends TipPartialsModel {

    public static final String SEQUENCE_ERROR_MODEL = "sequenceErrorModel";
    public static final String BASE_ERROR_RATE = "baseErrorRate";
    public static final String AGE_RELATED_RATE = "ageRelatedErrorRate";
    public static final String ERROR_TYPE1 = "type1";
    public static final String ERROR_TYPE2 = "type2";

    public static final String EXCLUDE = "exclude";
    public static final String INCLUDE = "include";

    public SequenceErrorModel(TreeModel treeModel, TaxonList includeTaxa, TaxonList excludeTaxa,
                              Parameter baseErrorRateParameter, Parameter ageRelatedErrorRateParameter) {
        super(SEQUENCE_ERROR_MODEL, treeModel, includeTaxa, excludeTaxa);

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
        double base = 0.0;
        double rate = 0.0;

        if (baseErrorRateParameter != null) {
            base = baseErrorRateParameter.getParameterValue(0);
        }

        if (ageRelatedErrorRateParameter != null) {
            rate = ageRelatedErrorRateParameter.getParameterValue(0);
        }

        int[] states = this.states[nodeIndex];

        if (!excluded[nodeIndex]) {
            double pUndamaged = (1.0 - base);

            if (ageRelatedErrorRateParameter != null) {
                double age = treeModel.getNodeHeight(treeModel.getExternalNode(nodeIndex));
                pUndamaged *= Math.exp(-rate * age);
            }

            int k = 0;
            for (int j = 0; j < patternCount; j++) {
                switch (states[j]) {
                    case 0: // is an A
                        partials[k] = pUndamaged;
                        partials[k + 1] = 0.0;
                        partials[k + 2] = 1.0 - pUndamaged;
                        partials[k + 3] = 0.0;
                        break;
                    case 1: // is an C
                        partials[k] = 0.0;
                        partials[k + 1] = pUndamaged;
                        partials[k + 2] = 0.0;
                        partials[k + 3] = 1.0 - pUndamaged;
                        break;
                    case 2: // is an G
                        partials[k] = 1.0 - pUndamaged;
                        partials[k + 1] = 0.0;
                        partials[k + 2] = pUndamaged;
                        partials[k + 3] = 0.0;
                        break;
                    case 3: // is an T
                        partials[k] = 0.0;
                        partials[k + 1] = 1.0 - pUndamaged;
                        partials[k + 2] = 0.0;
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

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public String getParserName() { return SEQUENCE_ERROR_MODEL; }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            TreeModel treeModel = (TreeModel)xo.getChild(TreeModel.class);

            Parameter baseDamageRateParameter = null;
            if (xo.hasChildNamed(BASE_ERROR_RATE)) {
                baseDamageRateParameter = (Parameter)xo.getElementFirstChild(BASE_ERROR_RATE);
            }

            Parameter ageRelatedRateParameter = null;
            if (xo.hasChildNamed(AGE_RELATED_RATE)) {
                ageRelatedRateParameter = (Parameter)xo.getElementFirstChild(AGE_RELATED_RATE);
            }

            if (baseDamageRateParameter == null && ageRelatedRateParameter == null) {
                throw new XMLParseException("You must specify one or other or both of " + BASE_ERROR_RATE + " and " + AGE_RELATED_RATE + " parameters");
            }

            TaxonList includeTaxa = null;
            TaxonList excludeTaxa = null;

            if (xo.hasChildNamed(INCLUDE)) {
                includeTaxa = (TaxonList)xo.getElementFirstChild(INCLUDE);
            }

            if (xo.hasChildNamed(EXCLUDE)) {
                excludeTaxa = (TaxonList)xo.getElementFirstChild(EXCLUDE);
            }

            SequenceErrorModel aDNADamageModel =  new SequenceErrorModel(treeModel, includeTaxa, excludeTaxa, baseDamageRateParameter, ageRelatedRateParameter);

            System.out.println("Using aDNA damage model.");

            return aDNADamageModel;
        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public String getParserDescription() {
            return
                    "This element returns a model that allows for post-mortem DNA damage.";
        }

        public Class getReturnType() { return SequenceErrorModel.class; }

        public XMLSyntaxRule[] getSyntaxRules() { return rules; }

        private XMLSyntaxRule[] rules = new XMLSyntaxRule[] {
                new ElementRule(TreeModel.class),
                new ElementRule(BASE_ERROR_RATE, Parameter.class, "The base error rate per site per sequence", true),
                new ElementRule(AGE_RELATED_RATE, Parameter.class, "The error rate per site per unit time", true),
                new XORRule(
                        new ElementRule(INCLUDE, TaxonList.class, "A set of taxa to which to apply the damage model to"),
                        new ElementRule(EXCLUDE, TaxonList.class, "A set of taxa to which to not apply the damage model to")
                        , true)
        };
    };

    private final Parameter baseErrorRateParameter;
    private final Parameter ageRelatedErrorRateParameter;
}