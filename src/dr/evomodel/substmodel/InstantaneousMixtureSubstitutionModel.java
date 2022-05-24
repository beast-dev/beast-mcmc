package dr.evomodel.substmodel;

import dr.evolution.datatype.DataType;
import dr.evomodel.substmodel.codon.GY94CodonModel;
import dr.inference.glm.GeneralizedLinearModel;
import dr.inference.loggers.LogColumn;
import dr.inference.model.BayesianStochasticSearchVariableSelection;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.util.Citation;
import dr.util.CommonCitations;

import java.util.Collections;
import java.util.List;

public class InstantaneousMixtureSubstitutionModel extends ComplexSubstitutionModel implements ParameterReplaceableSubstitutionModel {


    public InstantaneousMixtureSubstitutionModel(String name,
                                DataType dataType,
                                FrequencyModel rootFreqModel,
                                SubstitutionRateMatrixMixture substitutionRateMatrixMixture,
                                Parameter mixtureWeights) {
//                                boolean normalizeWeights) {

        super("instantaneousMixtureSubstitutionModel", dataType, rootFreqModel, null);
        this.substitutionRateMatrixMixture = substitutionRateMatrixMixture;
        this.mixtureWeights = mixtureWeights;
        // TODO generalize
        this.numComponents = 2;
//        System.err.println("constructing InstantaneousMixtureSubstitutionModel");
        addVariable(mixtureWeights);
    }

    public double[] getRates() {
        double[][] rateComponents = substitutionRateMatrixMixture.getParameterAsMatrix();
        double p = mixtureWeights.getParameterValue(0);
        double[] w = new double[]{p, 1.0 - p};
//        System.err.println("Mixing with p = " + p);

        double[] rates = new double[rateComponents.length];
        for (int i = 0; i < rateComponents.length; i++) {
            for (int j = 0; j < numComponents; j++) {
                rates[i] += w[j] * rateComponents[i][j];
            }
            rates[i] = Math.exp(rates[i]);
        }
        return rates;
    }

    protected void setupRelativeRates(double[] rates) {
        System.arraycopy(getRates(),0,rates,0,rates.length);
    }

    @Override
    public String getDescription() {
        return "Substitution model from log-linear combinations of Q matrices.";
    }

    @Override
    public ParameterReplaceableSubstitutionModel factory(List<Parameter> oldParameters, List<Parameter> newParameters) {
        Parameter weights = mixtureWeights;
        for (int i = 0; i < oldParameters.size(); i++) {
            Parameter oldParameter = oldParameters.get(i);
            Parameter newParameter = newParameters.get(i);
            if (oldParameter == mixtureWeights) {
                weights = newParameter;
            } else {
                throw new RuntimeException("Parameter not found in InstantaneousMixtureSubstitutionModel.");
            }
        }
        return new InstantaneousMixtureSubstitutionModel(super.getModelName(), super.dataType, super.freqModel, substitutionRateMatrixMixture, weights);
    }

    private int numComponents;
    private SubstitutionRateMatrixMixture substitutionRateMatrixMixture;
    private Parameter mixtureWeights;
}
