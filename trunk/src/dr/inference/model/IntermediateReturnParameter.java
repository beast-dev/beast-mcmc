package dr.inference.model;

/**
 * Created by max on 10/22/14.
 */
//Designed to return a data matrix post computation if asked. Designed for latent liabilities
public class IntermediateReturnParameter extends CompoundParameter {
    IntermediateReturnLikelihood IRL;

    public IntermediateReturnParameter(IntermediateReturnLikelihood IRL) {
        super(null);
        this.IRL = IRL;
    }

    @Override
    public int getDimension(){
        return IRL.returnIntermediate().getDimension();
    }

    public Parameter getParameter(int PID) {
        return IRL.returnIntermediate(PID);
    }
};

