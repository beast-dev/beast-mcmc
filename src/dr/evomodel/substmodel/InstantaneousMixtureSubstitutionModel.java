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
                                Parameter mixtureWeights,
                                boolean transform) {

        super("instantaneousMixtureSubstitutionModel", dataType, rootFreqModel, null);

        this.substitutionModelList = substitutionModelList;
        this.mixtureWeights = mixtureWeights;
        this.alphabetSize = substitutionModelList.get(0).getFrequencyModel().getFrequencyCount();
        this.alphabetSize2 = alphabetSize * alphabetSize;
        this.nRates = alphabetSize * (alphabetSize - 1);
        this.numComponents = substitutionModelList.size();
        this.pOneMinusP = numComponents == 2 && mixtureWeights.getSize() == 1;
        this.transform = transform;

        if (!checkWeightDimension()) {throw new RuntimeException("Mismatch between number of mixture weights and number of substitution models.");}

        if (!checkStateSpaces(substitutionModelList)) {
            throw new RuntimeException("Not all substitution models in instantaneousMixtureSubstitutionModel have same state space size.");
        }

        addVariable(mixtureWeights);
        for (SubstitutionModel substitutionModel : substitutionModelList) {
            addModel(substitutionModel);
        }
    }

    private boolean checkWeightDimension() {
        boolean weightsValid = true;
        if (numComponents != mixtureWeights.getSize()) {
            if ((!pOneMinusP) || (transform && mixtureWeights.getSize() != numComponents - 1)) {
                weightsValid = false;
            }
        }
        return weightsValid;
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

    private double[] getMixtureProportions() {
        double[] w = mixtureWeights.getParameterValues();

        if ( transform ) {
            double[] tmp = new double[w.length + 1];
            tmp[0] = 1.0;
            double totalSum = 1.0;
            for (int i = 0; i < tmp.length - 1; i++) {
                tmp[i+1] = w[i];
                totalSum += w[i];
            }
            int donothing = 0;
            for (int i = 0; i < tmp.length; i++) {
                tmp[i] /= totalSum;
            }
            w = tmp;
        } else if ( pOneMinusP ) {
            double p = mixtureWeights.getParameterValue(0);
            w = new double[]{p, 1.0 - p};
        }

        if ( checkWeights ) {
            double s = 0.0;
            for ( int i = 0; i < numComponents; i++ ) {
                if ( w[i] < 0.0 || w[i] > 1.0 ) {
                    throw new RuntimeException("Mixing proportion " + i + " has value (" + w[i] + ") outside allowed range of [0,1]");
                }
                s += w[i];
            }
            if ( Math.abs(s - 1.0) > Double.MIN_VALUE ) {
                throw new RuntimeException("Mixing proportions do not sum to 1");
            }
        }

        return w;
    }

    public double[] getRates() {
        double[] rateComponents = getComponentRates();
        double[] w = getMixtureProportions();

        double[] rates = new double[nRates];
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
        return new InstantaneousMixtureSubstitutionModel(name, dt, frequencies, subsModels, weights, transform);
    }

    private int alphabetSize;
    private int alphabetSize2;
    private int nRates;
    private int numComponents;
    private List<SubstitutionModel> substitutionModelList;
    private boolean pOneMinusP;
    private boolean transform;
    private boolean checkWeights = true;
    private Parameter mixtureWeights;
}
