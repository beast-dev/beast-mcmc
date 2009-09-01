package dr.evomodel.substmodel;

import dr.evolution.datatype.Microsatellite;
import dr.inference.model.Parameter;

/**
 * @author Chieh-Hsi Wu
 *
 * An abstract for One Phase Microsatellite Models
 */
public abstract class OnePhaseModel extends MicrosatelliteModel{

    /**
     * Constructor
     *
     * @param name              Model name
     * @param microsatellite    Microsatellite data type
     * @param freqModel         Equilibrium frequencies
     * @param parameter         Infinitesimal rates
     */
    public OnePhaseModel(String name, Microsatellite microsatellite, FrequencyModel freqModel, Parameter parameter){
        super(name, microsatellite, freqModel, parameter);
    }

    /**
     * The One Phase Models are special cases of the birth-death chain,
     * and therefore we can use this to calculate the stationay distribution
     * given a infinitesimal rate matrix.
     */
    public void computeStationaryDistribution(){
        double[] pi = new double[stateCount];
        double[] q = getRates();
        int k = 0;
        pi[0] = 1.0;
        double piSum = 1.0;
        for(int i = 1; i < stateCount; i++){
            pi[i] = pi[i-1]*q[k]/q[k+stateCount-1];
            k = k+stateCount;
            piSum = piSum+pi[i];
        }

        for(int i = 0; i < stateCount; i++){
            pi[i] = pi[i]/piSum;

        }
        freqModel = new FrequencyModel(dataType,pi);        
        super.computeStationaryDistribution();
    }


}
