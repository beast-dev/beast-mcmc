package dr.evomodel.substmodel;


import dr.evolution.datatype.Microsatellite;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;

import java.util.ArrayList;

/**
 * @author Chieh-Hsi Wu
 *
 * Implements Two Phase Model of microsatellites
 */
public class TwoPhaseModel extends MicrosatelliteModel{


    /*
     *Parameters for setting up the infinitesimal rate matrix.
     */
    private Parameter geoParam;
    private Parameter onePhasePrParam;
    private Parameter transformParam;
    private boolean estimateSubmodelParams = false;
    private boolean useEmpiricalFreqs = false;
    private ArrayList<Parameter> submodelParameters = null;
    private boolean updateSubmodelRates = false;

    public static final String TWO_PHASE_MODEL = "TWOPHASEModel";

    /**
     * Constructor
     *
     * @param microsatellite            microsatellite data type
     * @param freqModel                 user specified initial equilibrium frequencies
     * @param submodel                  this submodel of this
     * @param onePhasePrParam           the probability that the rate will conform to the one phase model
     * @param geoParam                  the probability of a success in geometric distribution of the TwoPhaseModel
     * @param transParam                the parameter for transforming onePhasePrParam and geoParam when onePhasePrParam > 0
     * @param estimateSubmodelParams    indicate whether the submodel parameters are estimated
     */

    public TwoPhaseModel(
            Microsatellite microsatellite,
            FrequencyModel freqModel,
            OnePhaseModel submodel,
            Parameter onePhasePrParam,
            Parameter geoParam,
            Parameter transParam,
            boolean estimateSubmodelParams){

        super(TWO_PHASE_MODEL, microsatellite, freqModel, null);

        this.subModel = submodel;
        this.estimateSubmodelParams = estimateSubmodelParams;

        if(this.estimateSubmodelParams){
            submodelParameters = new ArrayList<Parameter>();
            for(int i = 0; i < subModel.getNestedParameterCount(); i++){
                addVariable(subModel.getNestedParameter(i));
                submodelParameters.add(subModel.getNestedParameter(i));
            }
            updateSubmodelRates = true;
        }

        this.geoParam = geoParam;
        this.onePhasePrParam = onePhasePrParam;

        addVariable(this.geoParam);
        addVariable(this.onePhasePrParam);

        this.estimateSubmodelParams = estimateSubmodelParams;

        if(transParam != null){
            this.transformParam = transParam;
            System.out.println("TwoPhaseModel transformParameter value: "+transformParam);
        }else{
            this.transformParam = new Parameter.Default(0.0);
            System.out.println("TwoPhaseModel - default transformParameter value : "+transformParam);
        }

        setupInfinitesimalRates();

        if(freqModel == null){
            useEmpiricalFreqs = true;
            System.out.println("TwoPhaseModel: use empirical frequencies");
        }else{
            useEmpiricalFreqs = false;
        }


    }

    private Parameter transOnePhase;
    private Parameter transGeo;
    private void transform(){
        double e = transformParam.getParameterValue(0);
        double m = geoParam.getParameterValue(0);
        double p = onePhasePrParam.getParameterValue(0);
        if(p < 1 - e && m < 1 - e || p ==m){
            transOnePhase = onePhasePrParam;
            transGeo = geoParam;
        }else if(m > Math.max(1 - e,p)){
            p = p*(m-(m+e-1)/e)/m+(m+e-1)/e;
            transOnePhase = new Parameter.Default(p);
        }else if(p > Math.max(1 - e,m)){
            m = m*(p-(p+e-1)/e)/p+(p+e-1)/e;
            transGeo = new Parameter.Default(m);
        }
    }


    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        if(submodelParameters !=null && submodelParameters.indexOf((Parameter)variable) != -1){
            updateSubmodelRates = true;
        }
        updateMatrix = true;

    }

    public void setupInfinitesimalRates(){
        if(updateSubmodelRates){
            subModel.setupInfinitesimalRates();
            updateSubmodelRates = false;
        }
        transform();

        double geoParameter = transGeo.getParameterValue(0);
        double p = transOnePhase.getParameterValue(0);

        double[] condProbNum = new double[stateCount];
        for(int i = 1; i < stateCount; i++){
            condProbNum[i] = geoParameter*Math.pow((1.0 - geoParameter),i-1);
        }
        double condGeo = 0.0;

        for(int i = 0; i < stateCount; i++){
            double expansionGeoDenom = 1-Math.pow((1.0-geoParameter),stateCount - 1 - i);
            double contractionGeoDenom = 1 - Math.pow((1.0-geoParameter), i);
            double rowSum = 0.0;
            double submodelRate = 0.0;

            for(int j = 0; j < stateCount; j++){
                if(j < i){
                    condGeo = condProbNum[Math.abs(i-j)]/contractionGeoDenom;
                    submodelRate = subModel.getRate(i,i-1);
                }else if(j > i){
                    submodelRate = subModel.getRate(i,i+1);
                    condGeo = condProbNum[Math.abs(i-j)]/expansionGeoDenom;
                }

                if(i != j){

                    if(i == j + 1 || i == j - 1){
                        infinitesimalRateMatrix[i][j]= submodelRate*(p + (1 - p)*condGeo);
                    }else {
                        infinitesimalRateMatrix[i][j] = submodelRate*(1 - p)*condGeo;
                    }
                    rowSum = rowSum+infinitesimalRateMatrix[i][j];
                }


            }

            infinitesimalRateMatrix[i][i] = 0.0-rowSum;
        }
    }

    public void setupMatrix(){
        setupInfinitesimalRates();
        super.setupMatrix();
    }

    public void computeStationaryDistribution() {

        if(useEmpiricalFreqs){
            setupEmpiricalStationaryFrequencies();
        }
        super.computeStationaryDistribution();
    }

    public MicrosatelliteModel getSubModel(){
        return subModel;
    }

    public Parameter getGeometricParamter(){
        return geoParam;
    }

    public Parameter getOnePhasePrParamter(){
        return onePhasePrParam;

    }

    public Parameter getTransGeometricParamter(){
        return transGeo;
    }

    public Parameter getTransOnePhasePrParamter(){
        return transOnePhase;
    }

    public Parameter getTransformParam(){
        return transformParam;
    }

    public boolean isEstimatingSubmodelParams(){
        return estimateSubmodelParams;
    }

}
