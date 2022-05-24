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
                                List<SubstitutionModel> substitutionModelList,
                                Parameter mixtureWeights) {
//                                boolean normalizeWeights) {

        super("instantaneousMixtureSubstitutionModel", dataType, rootFreqModel, null);
        this.substitutionModelList = substitutionModelList;
        this.mixtureWeights = mixtureWeights;
        // TODO generalize
        this.numComponents = 2;
//        System.err.println("constructing InstantaneousMixtureSubstitutionModel");
        addVariable(mixtureWeights);
        for (int i = 0; i < substitutionModelList.size(); i++) {
            addModel(substitutionModelList.get(i));
        }
    }

    private double[][] getRatesAsMatrix() {

        int alphabetSize = substitutionModelList.get(0).getFrequencyModel().getFrequencyCount();
        int fullSize = alphabetSize * alphabetSize;
        int size = alphabetSize * (alphabetSize - 1);
        int halfSize = size/2;

        double[][] instantaneousRateMatrices = new double[size][substitutionModelList.size()];

        for (int m = 0; m < substitutionModelList.size(); ++m) {
            // Get Q-matrix
            SubstitutionModel model = substitutionModelList.get(m);
            // Full matrix stored row-wise
            double[] qMatrix = new double[fullSize];
            model.getInfinitesimalMatrix(qMatrix);

            // get base freqs
            double[] freqs = model.getFrequencyModel().getFrequencies();
            for (int i = 0; i < alphabetSize; i++) {
                freqs[i] = Math.log(freqs[i]);
            }

            double[] rates = new double[size];
            int idx = 0;
            for (int i = 0; i < (alphabetSize - 1); i++) {
                for (int j = (i + 1); j < alphabetSize; j++) {
                    rates[idx] = Math.log(qMatrix[i*alphabetSize + j]) - freqs[j];
                    rates[idx + halfSize] = Math.log(qMatrix[j*alphabetSize + i]) - freqs[i];
                    idx++;
                }
            }
            instantaneousRateMatrices[m] = rates;
        }
        return instantaneousRateMatrices;
    }

    public double[] getRates() {
        double[][] rateComponents = getRatesAsMatrix();
        double p = mixtureWeights.getParameterValue(0);
        double[] w = new double[]{p, 1.0 - p};
//        System.err.println("Mixing with p = " + p);
//        System.err.println("rateComponents[0][0] = " + rateComponents[0][0] + "; rateComponents[0][1] = " + rateComponents[0][1]);
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
        return new InstantaneousMixtureSubstitutionModel(super.getModelName(), super.dataType, super.freqModel, substitutionModelList, weights);
    }

    private int numComponents;
    private List<SubstitutionModel> substitutionModelList;
    private Parameter mixtureWeights;
}
