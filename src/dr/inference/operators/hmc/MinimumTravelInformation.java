package dr.inference.operators.hmc;

public class MinimumTravelInformation {

    final double time;
    final int[] index;
    final AbstractParticleOperator.Type type;

    MinimumTravelInformation(double minTime, int[] minIndex, AbstractParticleOperator.Type type) { //todo: merge multiple constructors
        this.time = minTime;
        this.index = minIndex;
        this.type = type;
    }

    MinimumTravelInformation(double minTime, int[] minIndex) {
        this.time = minTime;
        this.index = minIndex;
        this.type = AbstractParticleOperator.Type.NONE;
    }

    public MinimumTravelInformation(double minTime, int minIndex) {
        this.time = minTime;
        this.index = new int[]{minIndex};
        this.type = AbstractParticleOperator.Type.NONE;
    }

    MinimumTravelInformation(double minTime, int minIndex, AbstractParticleOperator.Type type) {
        this.time = minTime;
        this.index = new int[]{minIndex};
        this.type = type;
    }

    @SuppressWarnings("unused")
    public MinimumTravelInformation(double minTime, int[] minIndex, int ordinal) {
        this(minTime, minIndex, AbstractParticleOperator.Type.castFromInt(ordinal));
    }

    public boolean equals(Object obj) {

        if (obj instanceof MinimumTravelInformation) {
            MinimumTravelInformation rhs = (MinimumTravelInformation) obj;
//            return Math.abs(this.time - rhs.time) < 1E-5 && this.index == rhs.index;
            return this.time == rhs.time && this.index == rhs.index;
        }
        return false;
    }

    public String toString() {
        return "time = " + time + " @ " + index;
    }
}
