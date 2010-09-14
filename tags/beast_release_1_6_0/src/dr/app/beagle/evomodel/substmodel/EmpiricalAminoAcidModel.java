package dr.app.beagle.evomodel.substmodel;

import dr.evolution.datatype.AminoAcids;
import dr.evomodel.substmodel.EmpiricalRateMatrix;
import dr.inference.model.Parameter;

/**
 * @author Marc A. Suchard
 * @author Alexei Drummond
 * @author Andrew Rambaut
 */
public class EmpiricalAminoAcidModel extends BaseSubstitutionModel {

    public EmpiricalAminoAcidModel(EmpiricalRateMatrix rateMatrix, FrequencyModel freqModel) {

        super(rateMatrix.getName(), rateMatrix.getDataType(), freqModel, null);

        if (freqModel == null) {
            double[] freqs = rateMatrix.getEmpiricalFrequencies();
            this.freqModel = new FrequencyModel(AminoAcids.INSTANCE, new Parameter.Default(freqs));
        }

        this.rateMatrix = rateMatrix;
    }

    protected void frequenciesChanged() {
        // Nothing to precalculate
    }

    protected void ratesChanged() {
        // Nothing to precalculate
    }

    protected void setupRelativeRates(double[] rates) {
        double[] empiricalRates = rateMatrix.getEmpiricalRates();
        System.arraycopy(empiricalRates, 0, rates, 0, empiricalRates.length);

    }

    private EmpiricalRateMatrix rateMatrix;
}
