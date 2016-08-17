package test.dr.evomodel.substmodel;

import test.dr.math.MathTestCase;
import dr.evomodel.substmodel.FrequencyModel;
import dr.evomodel.substmodel.ProductChainFrequencyModel;
import dr.evolution.datatype.Nucleotides;
import dr.math.matrixAlgebra.Vector;

import java.util.List;
import java.util.ArrayList;

/**
 * @author Marc A. Suchard
 */
public class ProductChainFrequencyModelTest extends MathTestCase {

    public void testFrequencyModel() {

        FrequencyModel firstPosition = new FrequencyModel(Nucleotides.INSTANCE, freq1);
        FrequencyModel secondPosition = new FrequencyModel(Nucleotides.INSTANCE, freq2);
        FrequencyModel thirdPosition = new FrequencyModel(Nucleotides.INSTANCE, freq3);

        List<FrequencyModel> freqModels = new ArrayList<FrequencyModel>(3);
        freqModels.add(firstPosition);
        freqModels.add(secondPosition);
        freqModels.add(thirdPosition);

        ProductChainFrequencyModel pcFreqModel = new ProductChainFrequencyModel("freq", freqModels);

        double[] freqs = pcFreqModel.getFrequencies();

        System.out.println("Freq.length = " + freqs.length);

        int pos1 = 2;
        int pos2 = 1;
        int pos3 = 3;

        int index = computeIndex(pos1, pos2, pos3);
        
        System.out.println("Entry: " + new Vector(pcFreqModel.decomposeEntry(index)));
        System.out.println("Freq = " + freqs[index]);
        System.out.println("Freq = " + computeFreq(pos1, pos2, pos3));

        assertEquals(computeFreq(pos1, pos2, pos3), freqs[index]);       
    }

    private int computeIndex(int i, int j, int k) {
        return i * 16 + j * 4 + k;
    }

    private double computeFreq(int i, int j, int k) {
        return freq1[i] * freq2[j] * freq3[k];
    }

    private static double[] freq1 = {0.3, 0.25, 0.20, 0.25}; // A,C,G,T
    private static double[] freq2 = {0.1, 0.4,  0.25, 0.25}; // A,C,G,T
    private static double[] freq3 = {0.3, 0.15, 0.25, 0.30}; // A,C,G,T
}
