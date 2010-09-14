package dr.evomodel.substmodel;

import dr.evolution.datatype.Microsatellite;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;

import java.util.ArrayList;


/**
 * @author Chieh-Hsi Wu
 *
 * An abstract for One Phase Microsatellite Models
 */
public abstract class OnePhaseModel extends MicrosatelliteModel{
    protected ArrayList<Variable<Double>> nestedParams = null;

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
        nestedParams=new ArrayList<Variable<Double>>();
    }

    /*
     * adding the parameters only if its not a submodel.
     */
    protected void addParam(Variable<Double> param){
        if(isNested){
            nestedParams.add(param);
        }else{
            super.addVariable(param);
        }
    }

    /*
     * get the parameters in this submodel
     */
    public Variable<Double> getNestedParameter(int i){
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
        if(useStationaryFreqs){
            computeOnePhaseStationaryDistribution();
        }
        super.computeStationaryDistribution();

    }


}
