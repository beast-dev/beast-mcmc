package dr.inference.operators.hmc;

public class CategoryBounceInformation {
    final double time;
    final int[] index;
    final AbstractParticleOperator.Type type;

    CategoryBounceInformation(double bounceTime, int[] bounceIndex) {
        this.time = bounceTime;
        this.index = bounceIndex;
        this.type = AbstractParticleOperator.Type.NONE;
    }
}
