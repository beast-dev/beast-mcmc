package dr.evomodel.substmodel;

import dr.inference.model.Parameter;
import dr.inference.model.Variable;
import dr.evolution.datatype.Microsatellite;

import java.util.Map;
import java.util.ArrayList;

/**
 * @author Chieh-Hsi Wu
 *
 * This class is used for model averaging over a subset of microsatellite models.
 */
public class MsatAveragingSubsetModel extends MicrosatelliteModel{
    private boolean logit;
    public static final int PROP_INDEX = 0;
    public static final int BIAS_CONST_INDEX = 1;
    public static final int BIAS_LIN_INDEX = 2;
    public static final int GEO_INDEX = 3;
    public static final double DEFAULT_PROP = 0.0;
    public static final double DEFAULT_BIAS_CONSTANT = 0.0;
    public static final double DEFAULT_BIAS_LINEAR = 0.0;
    public static final double DEFAULT_GEO = 1.0;
    public static final int PARAMETER_PRESENT = 1;

    public Parameter[][] paramModelMap;
    public Map<Integer, Integer> modelMap;
    public Parameter modelChoose;
    public Parameter modelIndicator;

    public ArrayList<Parameter> propRates = new ArrayList<Parameter>();
    public ArrayList<Parameter> biasConsts = new ArrayList<Parameter>();
    public ArrayList<Parameter> biasLins = new ArrayList<Parameter>();
    public ArrayList<Parameter> geos = new ArrayList<Parameter>();

    public MsatAveragingSubsetModel(
            Microsatellite msat,
            boolean logit,
            ArrayList<Parameter> propRates,
            ArrayList<Parameter> biasConsts,
            ArrayList<Parameter> biasLins,
            ArrayList<Parameter> geos,
            Parameter[][] paramModelMap,
            Parameter modelChoose,
            Parameter modelIndicator,
            Map<Integer, Integer> modelMap){

        super("MsatAveragingModel",msat, null, null);

        for(int i = 0; i < propRates.size(); i++){
            addVariable(propRates.get(i));
        }

        for(int i = 0; i < biasConsts.size(); i++){
            addVariable(biasConsts.get(i));

        }

        for(int i = 0; i < biasLins.size(); i++){
            addVariable(biasLins.get(i));
        }

        for(int i = 0; i < geos.size(); i++){
            addVariable(geos.get(i));
        }
        addVariable(modelChoose);
        addVariable(modelIndicator);
        this.propRates = propRates;
        this.biasConsts = biasConsts;
        this.biasLins = biasLins;
        this.geos = geos;

        this.logit = logit;
        this.modelChoose = modelChoose;
        this.modelIndicator = modelIndicator;
        this.paramModelMap = paramModelMap;
        this.modelMap = modelMap;

        setupInfinitesimalRates();
        setupMatrix();


    }

 

    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        modelUpdate = false;
        int paramIndex = -1;
        if(propRates.contains(variable)){
            paramIndex = PROP_INDEX;
        }else if(biasConsts.contains(variable)){
            paramIndex = BIAS_CONST_INDEX;
        }else if(biasLins.contains(variable)){
            paramIndex = BIAS_LIN_INDEX;
        }else if(geos.contains(variable)){
            paramIndex = GEO_INDEX;
        }
        //System.out.println("changeModel: "+ variable);
        if(paramIndex > -1){
            //System.out.println(modelMap);
            //System.out.println(getBitVectorValue());
            if(paramModelMap[paramIndex][modelMap.get(getBitVectorValue())] != null){
                updateMatrix = true;
                modelUpdate = true;

            }

        }else if(variable == modelChoose){
            updateMatrix = true;
            indicateModel();
            modelUpdate = true;

        }


    }

    public void indicateModel(){
        modelIndicator.setParameterValueQuietly(0, modelMap.get(getBitVectorValue()));
    }

    public void setupInfinitesimalRates(){
        double rowSum;
        double prop = DEFAULT_PROP;
        double biasConst = DEFAULT_BIAS_CONSTANT;
        double biasLin = DEFAULT_BIAS_LINEAR;
        double geo = DEFAULT_GEO;
        infinitesimalRateMatrix = new double[stateCount][stateCount];

        //if the rate proportional parameter is present
        if((int)modelChoose.getParameterValue(PROP_INDEX) == PARAMETER_PRESENT){
            //
            prop = getModelParameterValue(PROP_INDEX);
        }
        for(int i = 0; i < stateCount;i++){
            rowSum = 0.0;
            if(i - 1 > -1){
                infinitesimalRateMatrix[i][i - 1] = 1+prop*i;
                rowSum = rowSum + infinitesimalRateMatrix[i][i - 1];
            }

            if(i + 1 < stateCount){
                infinitesimalRateMatrix[i][i + 1] = 1+prop*i;
                rowSum = rowSum + infinitesimalRateMatrix[i][i + 1];
            }
            infinitesimalRateMatrix[i][i] = -rowSum;

        }

        if((int)modelChoose.getParameterValue(BIAS_CONST_INDEX) == PARAMETER_PRESENT){

            biasConst = getModelParameterValue(BIAS_CONST_INDEX);
            //System.out.print("biasConst: "+biasConst);
            if((int)modelChoose.getParameterValue(BIAS_LIN_INDEX) == PARAMETER_PRESENT){
                biasLin = getModelParameterValue(BIAS_LIN_INDEX);
                //System.out.print("biasLin: "+biasLin);
            }

            double[][] subRates =  infinitesimalRateMatrix;
            infinitesimalRateMatrix = new double[stateCount][stateCount];

            LinearBiasModel.setupInfinitesimalRates(
                    infinitesimalRateMatrix,
                    subRates,
                    biasConst,
                    biasLin,
                    stateCount,
                    logit
            );
        }

        if((int)modelChoose.getParameterValue(GEO_INDEX) == PARAMETER_PRESENT){
            geo = getModelParameterValue(GEO_INDEX);
            double[][] subRates =  infinitesimalRateMatrix;
            infinitesimalRateMatrix = new double[stateCount][stateCount];

            TwoPhaseModel.setupInfinitesimalRates(
                stateCount,
                geo,
                0.0,
                infinitesimalRateMatrix,
                subRates
            );
        }
    }

    private double getModelParameterValue(int paramIndex){
        int modelCode = modelMap.get(getBitVectorValue());
        return (paramModelMap[paramIndex][modelCode]).getParameterValue(0);

    }

    private int getBitVectorValue(){
        String bitVec = "";
        for(int i = 0; i < modelChoose.getDimension(); i++){
            bitVec = bitVec + (int)modelChoose.getParameterValue(i);
        }
        int bitVecVal = Integer.parseInt(bitVec,2);

        return bitVecVal;
    }

    public void computeStationaryDistribution(){
        if((int)modelChoose.getParameterValue(GEO_INDEX) == PARAMETER_PRESENT){
            computeTwoPhaseStationaryDistribution();
        }else{
            computeOnePhaseStationaryDistribution();
        }
        super.computeStationaryDistribution();
    }

}
