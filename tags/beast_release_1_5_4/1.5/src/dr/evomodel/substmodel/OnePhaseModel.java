package dr.evomodel.substmodel;

import dr.evolution.datatype.Microsatellite;
import dr.inference.model.Parameter;

import java.util.ArrayList;


/**
 * @author Chieh-Hsi Wu
 *
 * An abstract for One Phase Microsatellite Models
 */
public abstract class OnePhaseModel extends MicrosatelliteModel{
    protected ArrayList<Parameter> nestedParams = null;

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
        nestedParams=new ArrayList<Parameter>();
    }

    /*
     * adding the parameters only if its not a submodel.
     */
    protected void addParam(Parameter param){
        if(isNested){
            nestedParams.add(param);
        }else{
            super.addVariable(param);
        }
    }

    /*
     * get the parameters in this submodel
     */
    public Parameter getNestedParameter(int i){
        return nestedParams.get(i);
    }

    /*
     * get number of nested parameters in the submodel
     */
    public int getNestedParameterCount(){
        return nestedParams.size();
    }

    /*
     * The One Phase Models are special cases of the birth-death chain,
     * and therefore we can use this to calculate the stationay distribution
     * given a infinitesimal rate matrix.
     */
    public void computeStationaryDistribution(){
        double[] pi = new double[stateCount];

        pi[0] = 1.0;
        double piSum = 1.0;
        for(int i = 1; i < stateCount; i++){
            pi[i] = pi[i-1]*infinitesimalRateMatrix[i-1][i]/infinitesimalRateMatrix[i][i-1];
            piSum = piSum+pi[i];
        }

        for(int i = 0; i < stateCount; i++){
            pi[i] = pi[i]/piSum;

        }
        freqModel = new FrequencyModel(dataType,pi);
        super.computeStationaryDistribution();

    }


}
