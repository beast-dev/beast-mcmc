package dr.inference.operators.hmc;

public class MinimumTravelInformationBinary {

    final double time;
    final int index;
    final AbstractParticleOperator.Type type;

    MinimumTravelInformationBinary(double minTime, int minIndex, AbstractParticleOperator.Type type) {
        this.time = minTime;
        this.index = minIndex;
        this.type = type;
    }

    MinimumTravelInformationBinary(double minTime, int minIndex) {
        this.time = minTime;
        this.index = minIndex;
        this.type = AbstractParticleOperator.Type.NONE;
    }

    @SuppressWarnings("unused")
    public MinimumTravelInformationBinary(double minTime, int minIndex, int ordinal) {
        this(minTime, minIndex, AbstractParticleOperator.Type.castFromInt(ordinal));
    }

    public boolean equals(Object obj) {

        if (obj instanceof MinimumTravelInformationBinary) {
            MinimumTravelInformationBinary rhs = (MinimumTravelInformationBinary) obj;
//            return Math.abs(this.time - rhs.time) < 1E-5 && this.index == rhs.index;
            return this.time == rhs.time && this.index == rhs.index;
        }
        return false;
    }

    public String toString() {
        return "time = " + time + " @ " + index;
    }
}
