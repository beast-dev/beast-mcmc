package dr.evomodel.continuous;

import dr.inference.model.Parameter;
import dr.xml.*;
/**
 * Created with IntelliJ IDEA.
 * User: max
 * Date: 10/28/13
 * Time: 1:58 PM
 * To change this template use File | Settings | File Templates.
 */
public class LatentFactorModel {
    private Parameter data;
    private Parameter factor;
    private Parameter latent;
    public LatentFactorModel(Parameter dataIn, Parameter factorIn, Parameter latentIn){
        data=dataIn;
        factor=factorIn;
        latent=latentIn;
    }
}
