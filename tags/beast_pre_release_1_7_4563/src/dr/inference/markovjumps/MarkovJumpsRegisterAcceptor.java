package dr.inference.markovjumps;

import dr.inference.model.Parameter;

/**
 * @author Marc A. Suchard
 */

public interface MarkovJumpsRegisterAcceptor {

    public void addRegister(Parameter addRegisterParameter,
                            MarkovJumpsType type,
                            boolean scaleByTime);
}
