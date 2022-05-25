package dr.evomodel.substmodel;

import dr.evolution.datatype.DataType;
import dr.evomodel.substmodel.codon.GY94CodonModel;
import dr.inference.glm.GeneralizedLinearModel;
import dr.inference.loggers.LogColumn;
import dr.inference.model.BayesianStochasticSearchVariableSelection;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.math.matrixAlgebra.Vector;
import dr.util.Citation;
import dr.util.CommonCitations;

import javax.xml.crypto.Data;
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
        this.alphabetSize = substitutionModelList.get(0).getFrequencyModel().getFrequencyCount();
        this.alphabetSize2 = alphabetSize * alphabetSize;
        this.nRates = alphabetSize * (alphabetSize - 1);
        // TODO generalize
        this.numComponents = 2;

        if (!checkStateSpaces(substitutionModelList)) {
            throw new RuntimeException("Not all substitution models in instantaneousMixtureSubstitutionModel have same state space size.");
        }

        addVariable(mixtureWeights);
        for (SubstitutionModel substitutionModel : substitutionModelList) {
            addModel(substitutionModel);
        }
    }

    private boolean checkStateSpaces(List<SubstitutionModel> substitutionModelList) {
        boolean sameSize = true;
        for (SubstitutionModel substitutionModel : substitutionModelList) {
            if ( substitutionModel.getFrequencyModel().getFrequencyCount() != alphabetSize) {
                sameSize = false;
            }
        }
        return sameSize;
    }

    private double[] getComponentRates() {

        int halfSize = nRates/2;

        // 1:nRates are model 0, size+(1:nRates) are model 1, and so on
        double[] rates = new double[numComponents*nRates];
        for (int m = 0; m < numComponents; ++m) {
            // Get Q-matrix
            SubstitutionModel model = substitutionModelList.get(m);
            // Full matrix stored row-wise
            double[] qMatrix = new double[alphabetSize2];
            model.getInfinitesimalMatrix(qMatrix);

            // get base freqs
            double[] freqs = model.getFrequencyModel().getFrequencies();
            for (int i = 0; i < alphabetSize; i++) {
                freqs[i] = Math.log(freqs[i]);
            }

            int idx = 0;
            for (int i = 0; i < (alphabetSize - 1); i++) {
                for (int j = (i + 1); j < alphabetSize; j++) {
                    rates[idx + m * nRates] = Math.log(qMatrix[i*alphabetSize + j]) - freqs[j];
                    rates[idx + halfSize + m * nRates] = Math.log(qMatrix[j*alphabetSize + i]) - freqs[i];
                    idx++;
                }
            }
        }
        return rates;
    }

    public double[] getRates() {
        double[] rateComponents = getComponentRates();

        double p = mixtureWeights.getParameterValue(0);
        double[] w = new double[]{p, 1.0 - p};
//        System.err.println("Mixing with p = " + p);
        if (p < 0.0 || p > 1.0) {
            throw new RuntimeException("Mixing proportion is " + p + " outside allowed range of [0,1]");
        }
//        System.err.println("rateComponents[0][0] = " + rateComponents[0][0] + "; rateComponents[0][1] = " + rateComponents[0][1]);
        double[] rates = new double[nRates];
        int fuck = 0;
        for (int i = 0; i < nRates; i++) {
            for (int j = 0; j < numComponents; j++) {
                rates[i] += w[j] * rateComponents[i + j * nRates];
            }
            rates[i] = Math.exp(rates[i]);
        }
//        System.err.println(new Vector(rates));
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
        List<SubstitutionModel> subsModels = substitutionModelList;
        FrequencyModel frequencies = freqModel;
        DataType dt = dataType;
        String name = super.getModelName();

        for (int i = 0; i < oldParameters.size(); i++) {
            Parameter oldParameter = oldParameters.get(i);
            Parameter newParameter = newParameters.get(i);
            if (oldParameter == mixtureWeights) {
                weights = newParameter;
            } else {
                throw new RuntimeException("Parameter not found in InstantaneousMixtureSubstitutionModel.");
            }
        }
        return new InstantaneousMixtureSubstitutionModel(name, dt, frequencies, subsModels, weights);
    }

    private int alphabetSize;
    private int alphabetSize2;
    private int nRates;
    private int numComponents;
    private List<SubstitutionModel> substitutionModelList;
    private Parameter mixtureWeights;
}
