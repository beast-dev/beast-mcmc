package dr.evomodel.treedatalikelihood.continuous.backprop;

import dr.evomodel.treedatalikelihood.continuous.ContinuousDataLikelihoodDelegate;
import dr.inference.model.Likelihood;
import dr.inference.model.Parameter;

public interface BackPropParameterProvider {

    /** The BEAST parameter that will receive a gradient. */
    Parameter getRawParameter();

    /** How to map primitive OU gradients into this parameter's entries. */
    SingleParameterOUPrimitiveGradientMapper getMapper();

    /** For wiring / consistency checks. */
    Likelihood getTreeDataLikelihood();
    ContinuousDataLikelihoodDelegate getContinuousDataLikelihoodDelegate();

    class Default implements BackPropParameterProvider {

        private final Parameter rawParameter;
        private final SingleParameterOUPrimitiveGradientMapper mapper;
        private final Likelihood treeDataLikelihood;
        private final ContinuousDataLikelihoodDelegate delegate;

        public Default(Likelihood treeDataLikelihood,
                                                ContinuousDataLikelihoodDelegate delegate,
                                                Parameter rawParameter,
                                                SingleParameterOUPrimitiveGradientMapper mapper) {

            if (rawParameter == null) {
                throw new IllegalArgumentException("rawParameter cannot be null");
            }
            if (mapper == null) {
                throw new IllegalArgumentException("mapper cannot be null");
            }
//            if (mapper.getParameter() != rawParameter) {
//                throw new IllegalArgumentException(
//                        "Mapper parameter and rawParameter must be the same instance");
//            } //TODO this check consistency; currently silenced since rawParameter can be a CompoundParameter version
            this.treeDataLikelihood = treeDataLikelihood;
            this.delegate = delegate;
            this.rawParameter = rawParameter;
            this.mapper = mapper;
        }

        @Override
        public Parameter getRawParameter() {
            return rawParameter;
        }

        @Override
        public SingleParameterOUPrimitiveGradientMapper getMapper() {
            return mapper;
        }

        @Override
        public Likelihood getTreeDataLikelihood() {
            return treeDataLikelihood;
        }

        @Override
        public ContinuousDataLikelihoodDelegate getContinuousDataLikelihoodDelegate() {
            return delegate;
        }
    }
}
