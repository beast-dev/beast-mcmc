package dr.inference.operators.hmc;

import dr.inference.hmc.GradientWrtParameterProvider;
import dr.inference.hmc.PrecisionColumnProvider;
import dr.inference.hmc.PrecisionMatrixVectorProductProvider;
import dr.inference.model.Parameter;
import dr.util.TaskPool;

/**
 * @author Zhenyu Zhang
 * @author Marc A. Suchard
 */
abstract class AbstractZigZagOperator extends AbstractParticleOperator {

    AbstractZigZagOperator(GradientWrtParameterProvider gradientProvider,
                           PrecisionMatrixVectorProductProvider multiplicationProvider,
                           PrecisionColumnProvider columnProvider,
                           double weight, Options runtimeOptions, Parameter mask,
                           int threadCount) {

        super(gradientProvider, multiplicationProvider, columnProvider, weight, runtimeOptions, mask);
        this.taskPool = (threadCount > 1) ? new TaskPool(gradientProvider.getDimension(), threadCount) : null;
    }

    final TaskPool taskPool;
}
