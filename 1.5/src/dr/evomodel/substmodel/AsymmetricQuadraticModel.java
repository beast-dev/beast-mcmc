package dr.evomodel.substmodel;


import dr.evolution.datatype.Microsatellite;
import dr.inference.model.Parameter;

/**
 * @author Chieh-Hsi Wu
 *
 * Implements the Asymmetric Quadratic Model
 */
public class AsymmetricQuadraticModel extends OnePhaseModel{

    public static final String ASYMQUAD_MODEL = "ASYMQUADModel";

    /*
     *Parameters for setting up the infinitesimal rate matrix.
     */
    private Parameter expanConst;
    private Parameter expanLin;
    private Parameter expanQuad;
    private Parameter contractConst;
    private Parameter contractLin;
    private Parameter contractQuad;


    /**
     * Constructor
     *
     * @param microsatellite    Mirosatellite data type
     * @param freqModel         Frequency model
     */
    public AsymmetricQuadraticModel(Microsatellite microsatellite, FrequencyModel freqModel){
        this(
                microsatellite,
                freqModel,
                null, null, null, null, null, null,
                false);

    }


    public AsymmetricQuadraticModel(Microsatellite microsatellite, FrequencyModel freqModel, boolean isNested){
        this(
                microsatellite,
                freqModel,
                null, null, null, null, null, null,
                isNested
        );

    }


    /**
     * Constructor
     *
     * @param microsatellite    Mirosatellite data type
     * @param freqModel         Frequency model
     * @param expanConst        Expansion constant
     * @param expanLin          Expansion linear coefficient
     * @param expanQuad         Expansion quadratic coefficient
     * @param contractConst     Contraction constant
     * @param contractLin       Contraction linear coefficient
     * @param contractQuad      Contraction quadratic coefficient
     * @param isNested          boolean indicating whether this object is a submodel of another microsatellite model
     */
    public AsymmetricQuadraticModel(Microsatellite microsatellite, FrequencyModel freqModel,
                    Parameter expanConst, Parameter expanLin, Parameter expanQuad,
                    Parameter contractConst, Parameter contractLin, Parameter contractQuad,
                    boolean isNested){

        super(ASYMQUAD_MODEL, microsatellite, freqModel,null);


        //The default setting of the parameters gives the same infinitesimal rates
        // as the StepwiseMutaionalModel class.
        this.expanConst = overrideDefault(new Parameter.Default(1.0), expanConst);
        this.expanLin = overrideDefault(new Parameter.Default(0.0), expanLin);
        this.expanQuad = overrideDefault(new Parameter.Default(0.0), expanQuad);
        this.contractConst = overrideDefault(this.expanConst, contractConst);
        this.contractLin = overrideDefault(this.expanLin, contractLin);
        this.contractQuad = overrideDefault(this.expanQuad, contractQuad);
        this.isNested = isNested;
        addParameters();

        printDetails();

        setupInfinitesimalRates();

        //calculate the default frequencies when not provieded by the user.
        if(freqModel == null){
            System.out.println("Creating AysmmetricQuadraticModel: using empirical frequencies");
            computeStationaryDistribution();
        }else{
            this.freqModel = freqModel;
        }

        addModel(this.freqModel);

    }

    private void addParameters(){
        addParam(this.expanConst);
        addParam(this.expanLin);
        addParam(this.expanQuad);
        if(this.contractConst != this.expanConst)
            addParam(this.contractConst);
        if(this.contractLin != this.expanLin)
            addParam(this.contractLin);
        if(this.contractQuad != this.expanQuad)
            addParam(this.contractQuad);
    }




    /*
     *  This method will override the default value of the parameter using the value specified by the user.
     */
    private Parameter overrideDefault(Parameter defaultParam, Parameter providedParam){
        if(providedParam != null && providedParam != defaultParam)
            return providedParam;
        return defaultParam;
    }


    /*
     * Setting up the infinitesimal Rates
     * The rates are defined by the following equations:
     * X -> X + 1 at rate u0 + u1(X - k) + u2(X - k)^2
     * X -> X - 1 at rate d0 + d1(X - k) + d2(X - k)^2
     */
    public void setupInfinitesimalRates(){


        double u0 = expanConst.getParameterValue(0);
        double u1 = expanLin.getParameterValue(0);
        double u2 = expanQuad.getParameterValue(0);

        double d0 = contractConst.getParameterValue(0);
        double d1 = contractLin.getParameterValue(0);
        double d2 = contractQuad.getParameterValue(0);

        double rowSum;
        for(int i = 0; i < stateCount;i++){
            rowSum = 0.0;
            if(i - 1 > -1){
                infinitesimalRateMatrix[i][i - 1] =d0+d1*i+d2*i*i;
                rowSum = rowSum + infinitesimalRateMatrix[i][i - 1];

            }

            if(i + 1 < stateCount){
                infinitesimalRateMatrix[i][i + 1] = u0+u1*i+u2*i*i;
                rowSum = rowSum + infinitesimalRateMatrix[i][i + 1];

            }

            infinitesimalRateMatrix[i][i] = rowSum*-1;

        }


    }


    public void setupMatrix(){
        setupInfinitesimalRates();
        super.setupMatrix();
    }

    public Parameter getExpansionConstant(){
        return expanConst;
    }

    public Parameter getExpansionLinear(){
        return expanLin;
    }

    public Parameter getExpansionQuad(){
        return expanQuad;
    }

    public Parameter getContractionConstant(){
        return contractConst;
    }

    public Parameter getContractionLinear(){
        return contractLin;
    }

    public Parameter getContractionQuad(){
        return contractQuad;
    }

    public void printDetails(){
        System.out.println("\n");
        System.out.println("Details of the asymmetric quadratic model and its parameters:");
        System.out.println("expansion constant:   "+expanConst.getParameterValue(0));
        System.out.println("expansion linear:     "+ expanLin.getParameterValue(0));
        System.out.println("expansion quadratic:  "+expanQuad.getParameterValue(0));
        System.out.println("contraction constant: "+contractConst.getParameterValue(0));
        System.out.println("contraction linear:   "+contractLin.getParameterValue(0));
        System.out.println("contraction quadratc: "+contractQuad.getParameterValue(0));
        System.out.println("a submodel:           "+isNested);
        System.out.println("\n");
    }

}
