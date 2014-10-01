package dr.evomodel.treelikelihood;

import dr.evolution.datatype.Nucleotides;
import dr.evolution.util.TaxonList;
import dr.evomodelxml.treelikelihood.SequenceErrorModelParser;
import dr.inference.model.Parameter;
import dr.inference.model.Statistic;

/**
 * This class incorporates uncertainty in the state at the tips of the tree and can
 * be used to model processes like sequencing error and DNA damage. It can have a fixed
 * (per site) base error rate and/or a time dependent error for which the probability
 * of no error decays over sampling time exponentially with a given rate. This model
 * is inspired by a brief description in Joe Felsenstein's book 'Inferring phylogenies'
 * (2004: Sinauer Associates) and was elaborated on for DNA damage in Rambaut et al
 * (2008, MBE
 * @author Andrew Rambaut
 * @version $Id$
 */
public class SequenceErrorModel extends TipStatesModel {
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
                              ErrorType errorType, Parameter baseErrorRateParameter,
                              Parameter ageRelatedErrorRateParameter,
                              Parameter indicatorParameter) {
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

        if (indicatorParameter != null) {
            this.indicatorParameter = indicatorParameter;
            addVariable(indicatorParameter);
        } else {
            this.indicatorParameter = null;
        }

        if (indicatorParameter != null) {
            addStatistic(new TaxonHasErrorsStatistic());
        }
    }

    protected void taxaChanged() {
        if (indicatorParameter != null && indicatorParameter.getDimension() <= 1) {
            this.indicatorParameter.setDimension(tree.getExternalNodeCount());
        }
    }

    @Override
    public Type getModelType() {
        return Type.PARTIALS;
    }

    @Override
    public void getTipStates(int nodeIndex, int[] tipStates) {
        throw new IllegalArgumentException("This model emits only tip partials");
    }

    @Override
    public void getTipPartials(int nodeIndex, double[] partials) {

        int[] states = this.states[nodeIndex];
        if (indicatorParameter == null || indicatorParameter.getParameterValue(nodeIndex) > 0.0) {

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
                    case Nucleotides.A_STATE: // is an A
                        partials[k] = pUndamaged;
                        partials[k + 1] = pDamagedTV;
                        partials[k + 2] = pDamagedTS;
                        partials[k + 3] = pDamagedTV;
                        break;
                    case Nucleotides.C_STATE: // is an C
                        partials[k] = pDamagedTV;
                        partials[k + 1] = pUndamaged;
                        partials[k + 2] = pDamagedTV;
                        partials[k + 3] = pDamagedTS;
                        break;
                    case Nucleotides.G_STATE: // is an G
                        partials[k] = pDamagedTS;
                        partials[k + 1] = pDamagedTV;
                        partials[k + 2] = pUndamaged;
                        partials[k + 3] = pDamagedTV;
                        break;
                    case Nucleotides.UT_STATE: // is an T
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
        } else {
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
    public class TaxonHasErrorsStatistic extends Statistic.Abstract {

        public TaxonHasErrorsStatistic() {
            super("hasErrors");
        }

        public int getDimension() {
            if (indicatorParameter == null) return 0;
            return indicatorParameter.getDimension();
        }

        public String getDimensionName(int dim) {
            return taxonMap.get(dim);
        }

        public double getStatisticValue(int dim) {
            return indicatorParameter.getParameterValue(dim);
        }

    }


    private final ErrorType errorType;
    private final Parameter baseErrorRateParameter;
    private final Parameter ageRelatedErrorRateParameter;
    private final Parameter indicatorParameter;
}