package dr.app.beagle.evomodel.substmodel;

import dr.inference.model.Parameter;

import java.util.List;

/**
 * @author Marc A. Suchard
 * @author Vladimir Minin
 *         <p/>
 *         A class for implementing a kronecker sum of CTMC models in BEAST using BEAGLE
 *         This work is supported by NSF grant 0856099
 *         <p/>
 *         O'Brien JD, Minin VN and Suchard MA (2009) Learning to count: robust estimates for labeled distances between
 *         molecular sequences. Molecular Biology and Evolution, 26, 801-814
 */

public class MarkovModulatedFrequencyModel extends FrequencyModel {

    public MarkovModulatedFrequencyModel(String name, List<FrequencyModel> freqModels) {
        super(name);
        this.freqModels = freqModels;
        int freqCount = 0;
        stateCount = freqModels.get(0).getFrequencyCount();
        numBaseModel = freqModels.size();
        for (int i = 0; i < numBaseModel; i++) {
            int size = freqModels.get(i).getFrequencyCount();
            if (stateCount != size) {
                throw new RuntimeException("MarkovModulatedFrequencyModel requires all frequencies model to have the same dimension");
            }           
            freqCount += size;
        }
        totalFreqCount = freqCount;
    }

    public void setFrequency(int i, double value) {
        throw new RuntimeException("Not implemented");
    }

    public double getFrequency(int index) {
        int whichModel = index / stateCount;
        int whichState = index % stateCount;
        
        // Assumes that all hidden states are equiprobable a priori
        return freqModels.get(whichModel).getFrequency(whichState) / numBaseModel;
    }

    public int getFrequencyCount() {
        return totalFreqCount;
    }

    public Parameter getFrequencyParameter() {
        throw new RuntimeException("Not implemented");
    }
     
    private List<FrequencyModel> freqModels;

    private final int numBaseModel;
    private final int totalFreqCount;
    private final int stateCount;
}
