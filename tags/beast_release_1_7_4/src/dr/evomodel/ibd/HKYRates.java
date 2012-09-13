package dr.evomodel.ibd;

import dr.evomodel.substmodel.FrequencyModel;
import dr.evomodel.substmodel.HKY;
import dr.inference.model.Parameter;

/**
 * Package: dr.evomodel.ibd
 * Description:
 * <p/>
 * <p/>
 * Created by
 * avaleks (alexander.alekseyenko@gmail.com)
 * Date: 05-Aug-2008
 * Time: 09:32:35
 */
public class HKYRates extends HKY {
    /**
     * Constructor
     */
    public HKYRates(Parameter kappaParameter, FrequencyModel freqModel) {
        super(kappaParameter, freqModel);
    }

    public double[] getRelativeRates(double[] rateMatrix) {
        double kappa = getKappa();
        double[] freq = getFrequencyModel().getFrequencies();
        // A - C - G - T
        rateMatrix[0] = -(freq[1] + freq[3]) - freq[2] * kappa;
        rateMatrix[1] = freq[1];
        rateMatrix[2] = freq[2] * kappa;
        rateMatrix[3] = freq[3];

        rateMatrix[4] = freq[0];
        rateMatrix[5] = -(freq[0] + freq[2]) - freq[3] * kappa;
        rateMatrix[6] = freq[2];
        rateMatrix[7] = freq[3] * kappa;

        rateMatrix[8] = freq[0] * kappa;
        rateMatrix[9] = freq[1];
        rateMatrix[10] = -(freq[1] + freq[3]) - freq[0] * kappa;
        rateMatrix[11] = freq[3];

        rateMatrix[12] = freq[0];
        rateMatrix[13] = freq[1] * kappa;
        rateMatrix[14] = freq[2];
        rateMatrix[15] = -(freq[0] + freq[1]) - freq[1] * kappa;


        return rateMatrix;
    }
}
