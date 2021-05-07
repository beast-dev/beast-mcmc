package dr.inference.operators.hmc;

public interface MassPreconditioningOptions {

    int preconditioningUpdateFrequency();
    int preconditioningDelay();
    int preconditioningMaxUpdate();
    int preconditioningMemory();

    public class Default implements MassPreconditioningOptions {
        final int preconditioningUpdateFrequency;
        final int preconditioningMaxUpdate;
        final int preconditioningDelay;
        final int preconditioningMemory;
        final boolean guessInitialMass;

        public Default(int preconditioningUpdateFrequency, int preconditioningMaxUpdate,
                       int preconditioningDelay, int preconditioningMemory, boolean guessInitialMass) {
            this.preconditioningUpdateFrequency = preconditioningUpdateFrequency;
            this.preconditioningMaxUpdate = preconditioningMaxUpdate;
            this.preconditioningDelay = preconditioningDelay;
            this.preconditioningMemory = preconditioningMemory;
            this.guessInitialMass = guessInitialMass;
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
    }
}
