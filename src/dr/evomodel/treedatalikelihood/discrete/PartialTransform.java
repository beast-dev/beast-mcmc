package dr.evomodel.treedatalikelihood.discrete;

public interface PartialTransform {

    String getName();

    void transform(double[] input, double[] output);

    PartialTransform IDENTITY = new PartialTransform() {
        public String getName() { return "identity"; }
        public void transform(double[] input, double[] output) {
            System.arraycopy(input, 0, output, 0, input.length);
        }
    };
}
