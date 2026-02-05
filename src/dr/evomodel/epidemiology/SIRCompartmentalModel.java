package dr.evomodel.epidemiology;

import dr.inference.model.Parameter;

public class SIRCompartmentalModel extends CompartmentalModel{

    Parameter transmissionRate;
    Parameter recoveryRate;
    Parameter samplingProportion;
    Parameter numS;
    Parameter numI;
    Parameter numR;
    Parameter origin;
    int numGridPoints;
    double cutOff;

    public SIRCompartmentalModel(
            Parameter transmissionRate,
            Parameter recoveryRate,
            Parameter samplingProportion,
            Parameter numS,
            Parameter numI,
            Parameter numR,
            Parameter origin,
            int numGridPoints,
            double cutOff){

        super("SIRCompartmentalModel");

        this.transmissionRate = transmissionRate;
        addVariable(transmissionRate);
        this.recoveryRate = recoveryRate;
        addVariable(recoveryRate);
        this.samplingProportion = samplingProportion;
        addVariable(samplingProportion);
        this.numS = numS;
        addVariable(numS);
        this.numI = numI;
        addVariable(numI);
        this.numR = numR;
        addVariable(numR);
        this.origin = origin;
        addVariable(origin);

        this.numGridPoints = numGridPoints;
        this.cutOff = cutOff;
    }

    public void simulateTrajectory(){

        // FILL IN



        // Code below is just for testing
        // for(int i = 0; i < numS.getDimension(); i++) {
        //     numS.setParameterValue(i, numS.getParameterValue(0) - i);
        // }
    }

}
