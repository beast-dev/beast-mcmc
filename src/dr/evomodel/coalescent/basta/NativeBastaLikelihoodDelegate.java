package dr.evomodel.coalescent.basta;

import dr.evomodel.treedatalikelihood.continuous.cdi.CDIJNIWrapper;
import dr.evomodel.treedatalikelihood.continuous.cdi.ResourceDetails;

/**
 * @author Marc A. Suchard
 */
public class NativeBastaLikelihoodDelegate extends BastaLikelihoodDelegate.AbstractBastaLikelihood {

    private final NativeBastaJniWrapper jni;

    public NativeBastaLikelihoodDelegate(String name) {
        super(name);
        jni = NativeBastaJniWrapper.getBastaJniWrapper();
    }

    @Override
    protected void functionOne() {
        jni.functionOne();
    }

    @Override
    protected void functionTwo() {
        jni.functionTwo();
    }
}
