package dr.inference.operators.hmc;

import dr.inference.model.Parameter;

public interface MassPreconditioningOptions {

    int preconditioningUpdateFrequency();
    int preconditioningDelay();
    int preconditioningMaxUpdate();
    int preconditioningMemory();
    Parameter preconditioningEigenLowerBound();
    Parameter preconditioningEigenUpperBound();

    class Default implements MassPreconditioningOptions {
        final int preconditioningUpdateFrequency;
        final int preconditioningMaxUpdate;
        final int preconditioningDelay;
        final int preconditioningMemory;
        final boolean guessInitialMass;
        final Parameter preconditioningEigenLowerBound;
        final Parameter preconditioningEigenUpperBound;

        public Default(int preconditioningUpdateFrequency, int preconditioningMaxUpdate,
                       int preconditioningDelay, int preconditioningMemory, boolean guessInitialMass,
                       Parameter eigenLowerBound, Parameter eigenUpperBound) {
            this.preconditioningUpdateFrequency = preconditioningUpdateFrequency;
            this.preconditioningMaxUpdate = preconditioningMaxUpdate;
            this.preconditioningDelay = preconditioningDelay;
            this.preconditioningMemory = preconditioningMemory;
            this.guessInitialMass = guessInitialMass;
            this.preconditioningEigenLowerBound = eigenLowerBound;
            this.preconditioningEigenUpperBound = eigenUpperBound;
        }

        @Override
        public int preconditioningUpdateFrequency() {
            return preconditioningUpdateFrequency;
        }

        @Override
        public int preconditioningDelay() {
            return preconditioningDelay;
        }

        @Override
        public int preconditioningMaxUpdate() {
            return preconditioningMaxUpdate;
        }

        @Override
        public int preconditioningMemory() {
            return preconditioningMemory;
        }

        @Override
        public Parameter preconditioningEigenLowerBound() {
            return preconditioningEigenLowerBound;
        }

        @Override
        public Parameter preconditioningEigenUpperBound() {
            return preconditioningEigenUpperBound;
        }
    }
}
