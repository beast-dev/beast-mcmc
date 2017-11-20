package dr.evomodel.treedatalikelihood.continuous.cdi;

/**
 * @author Marc A. Suchard
 */
public enum DriftType {

    NONE("no drift") {
        @Override
        public void doSomeWork(int dimTrait) {

        }

        @Override
        public boolean hasDrift() {
            return false;
        }
    },

    SIMPLE("simple drift") {
        @Override
        public void doSomeWork(int dimTrait) {

        }
    },

    SCALAR("proportional shear with drift") {
        @Override
        public void doSomeWork(int dimTrait) {

        }

        @Override
        public int getShearSize(int dim) {
            return 1;
        }
    },

    FULL("full shear matrix with drift") {
        @Override
        public void doSomeWork(int dimTrait) {

        }

        @Override
        public int getShearSize(int dim) {
            return dim * dim;
        }
    };

    private String name;

    DriftType(String name) {
        this.name = name;
    }

    abstract public void doSomeWork(int dimTrait);

    public int getShearSize(int dim) { return 0; }

    public boolean hasDrift() { return true; }

    public boolean hasShear(int dim) {
        return getShearSize(dim) != 0;
    }

    public String getName() { return name; }
}
