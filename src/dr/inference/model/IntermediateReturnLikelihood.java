package dr.inference.model;

/**
 * Created by max on 10/31/14.
 */
public abstract class IntermediateReturnLikelihood extends AbstractModelLikelihood{
    public IntermediateReturnLikelihood(String name){super(name);}

    Parameter returnIntermediate(int PID){return null;}
    Parameter returnIntermediate(){return null;}
}
