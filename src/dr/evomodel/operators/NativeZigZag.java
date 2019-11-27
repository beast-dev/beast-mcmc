package dr.evomodel.operators;

import dr.inference.hmc.PrecisionColumnProvider;
import dr.inference.model.Parameter;
import dr.inference.operators.hmc.AbstractParticleOperator;

public class NativeZigZag {

    private final PrecisionColumnProvider columnProvider;
    private final Parameter mask;

    public NativeZigZag(PrecisionColumnProvider columnProvider,
                        Parameter mask) {

        this.columnProvider = columnProvider;
        this.mask = mask;
    }

    private native void operate(double[] position,
                                double[] velocity,
                                double[] gradient,
                                double[] action,
                                double time);

    private native AbstractParticleOperator.MinimumTravelInformation[] getNextEvent(double[] position,
                                                                                    double[] velocity,
                                                                                    double[] gradient,
                                                                                    double[] action,
                                                                                    double time);
}
