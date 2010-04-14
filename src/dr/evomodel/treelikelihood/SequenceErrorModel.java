package dr.evomodel.treelikelihood;

import dr.evolution.util.TaxonList;
import dr.evomodelxml.treelikelihood.SequenceErrorModelParser;
import dr.inference.model.Parameter;

/**
 * This class incorporates uncertainty in the state at the tips of the tree and can
 * be used to model processes like sequencing error and DNA damage. It can have a fixed
 * (per site) base error rate and/or a time dependent error for which the probability
 * of no error decays over sampling time exponentially with a given rate. This model
 * is inspired by a brief description in Joe Felsenstein's book 'Inferring phylogenies'
 * (2004: Sinauer Associates) and was ellaborated on for DNA damage in Rambaut et al
 * (2008, MBE
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

    public SequenceErrorModel(TaxonList includeTaxa, TaxonList excludeTaxa,
                              ErrorType errorType, Parameter baseErrorRateParameter, Parameter ageRelatedErrorRateParameter) {
        super(SequenceErrorModelParser.SEQUENCE_ERROR_MODEL, includeTaxa, excludeTaxa);

        this.errorType = errorType;

        if (baseErrorRateParameter != null) {
            this.baseErrorRateParameter = baseErrorRateParameter;
            addVariable(this.baseErrorRateParameter);
        } else {
            this.baseErrorRateParameter = null;
        }

        if (ageRelatedErrorRateParameter != null) {
            this.ageRelatedErrorRateParameter = ageRelatedErrorRateParameter;
            addVariable(ageRelatedErrorRateParameter);
        } else {
            this.ageRelatedErrorRateParameter = null;
        }
    }

    protected void taxaChanged() {
        // nothing to do
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

    private final ErrorType errorType;
    private final Parameter baseErrorRateParameter;
    private final Parameter ageRelatedErrorRateParameter;
}