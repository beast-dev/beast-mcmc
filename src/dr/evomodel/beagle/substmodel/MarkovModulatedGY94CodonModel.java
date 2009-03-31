package dr.evomodel.beagle.substmodel;

import dr.evolution.datatype.HiddenCodons;
import dr.evolution.datatype.Codons;
import dr.inference.model.Parameter;
import dr.math.matrixAlgebra.Matrix;
import dr.math.matrixAlgebra.Vector;

/**
 * @author Marc A. Suchard
 */
public class MarkovModulatedGY94CodonModel extends GY94CodonModel {

    private static final byte RATE = 5;

    public MarkovModulatedGY94CodonModel(
            HiddenCodons codonDataType,
            Parameter omegaParameter,
            Parameter kappaParameter,
            Parameter switchingRates,
            FrequencyModel freqModel) {

        super(codonDataType, omegaParameter, kappaParameter, freqModel);

        this.hiddenClassCount = codonDataType.getHiddenClassCount();
        this.switchingRates = switchingRates;

    }

    protected void setupRelativeRates(double[] relativeRates) {
        double kappa = getKappa();
        double[] omega = omegaParameter.getParameterValues();
        double[] rates = switchingRates.getParameterValues();

        int stateCount = this.stateCount / hiddenClassCount;

        int index = 0;
        for (int i = 0; i < stateCount; i++) {
            for (int j = i + 1; j < stateCount; j++) {
                for (int h = 0; h < hiddenClassCount; h++) {

                    int d = getIndex(h * stateCount + i, h * stateCount + j, this.stateCount);
                    switch (rateMap[index]) {
                        case 0:
                            relativeRates[d] = 0.0;
                            break;                                // codon changes in more than one codon position
                        case 1:
                            relativeRates[d] = kappa;
                            break;                                // synonymous transition
                        case 2:
                            relativeRates[d] = 1.0;
                            break;                                // synonymous transversion
                        case 3:
                            relativeRates[d] = kappa * omega[h];
                            break;                                // non-synonymous transition
                        case 4:
                            relativeRates[d] = omega[h];
                            break;                                // non-synonymous transversion

                    }
                }
                index++;
            }
        }

//        System.err.println("rr = "+new Vector(relativeRates));
//        System.exit(-1);
        // Add the switching class rates
        double[] freqs = freqModel.getFrequencies();
        int rateIndex = 0;
        for (int g = 0; g < hiddenClassCount; g++) {
            for (int h = g + 1; h < hiddenClassCount; h++) {  // from g -> h
                for (int i = 0; i < stateCount; i++) {
                    int d = getIndex(g * stateCount + i, h * stateCount + i, this.stateCount);
                    // correct for the fact that setupMatrix post-multiplies these rates
//                    relativeRates[d] = rates[rateIndex] / freqs[i];
                }
                rateIndex++;
            }
        }
    }

//    public void setupRelativeRates() {
//        for (int i = 0; i < relativeRates.length; i++)
//            relativeRates[i] = 1.0;
//    }

    // Mapping: Matrix[i][j] = Compressed vector[i*(S - 3/2) - i^2 / 2 + j - 1]
    private static int getIndex(int i, int j, int S) {
        return (i * (2 * S - 3) - i * i) / 2 + j - 1;
    }
     
    protected void constructRateMap() {
        // Construct map for non-hidden states only
        hiddenClassCount = ((HiddenCodons) codonDataType).getHiddenClassCount();
        stateCount /= hiddenClassCount;
        super.constructRateMap();
        stateCount *= hiddenClassCount;
    }

    public static void main(String[] args) {
        GY94CodonModel codonModel = new GY94CodonModel(Codons.UNIVERSAL,
                new Parameter.Default(1.0), new Parameter.Default(1.0),
                new FrequencyModel(Codons.UNIVERSAL, new Parameter.Default(61,1.0/61.0)));
        EigenDecomposition ed1 = codonModel.getEigenDecomposition();
//        double[][] q = codonModel.getQ();
//        System.err.println("matrixQ = \n"+ new Matrix(q));

        MarkovModulatedGY94CodonModel mmCodonModel = new MarkovModulatedGY94CodonModel(HiddenCodons.UNIVERSAL_HIDDEN_2,
                new Parameter.Default(2,1.0), new Parameter.Default(2,1.0),
                new Parameter.Default(2.0),
                new FrequencyModel(HiddenCodons.UNIVERSAL_HIDDEN_2, new Parameter.Default(122,1.0/122.0)));
        EigenDecomposition ed2 = mmCodonModel.getEigenDecomposition();                

    }

    private int hiddenClassCount;
    private Parameter switchingRates;

}
