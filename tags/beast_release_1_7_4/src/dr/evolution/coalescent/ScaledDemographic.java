package dr.evolution.coalescent;

/**
 * Scale a demographic by a fixed factor.
 *
 * Minimal implementation, only value and intensity implemented.
 *
 * @author Joseph Heled
 *         Date: 19/11/2007
 */
public class ScaledDemographic extends DemographicFunction.Abstract {
    private final double scale;
    private final DemographicFunction demo;

    public ScaledDemographic(DemographicFunction demo, double scale) {
        super(demo.getUnits());
        this.scale = scale;
        this.demo = demo;
    }

    public double getDemographic(double t) {
        return demo.getDemographic(t) * scale;
    }

    public double getIntensity(double t) {
        return demo.getIntensity(t) / scale;
    }

    public double getIntegral(double start, double finish) {
        return (demo.getIntensity(finish) - demo.getIntensity(start)) / scale;
    }

    public double getInverseIntensity(double x) {
        throw new RuntimeException("unimplemented");
    }

    public int getNumArguments() {
        throw new RuntimeException("unimplemented");
    }

    public String getArgumentName(int n) {
        throw new RuntimeException("unimplemented");
    }

    public double getArgument(int n) {
        throw new RuntimeException("unimplemented");
    }

    public void setArgument(int n, double value) {
        throw new RuntimeException("unimplemented");
    }

    public double getLowerBound(int n) {
        throw new RuntimeException("unimplemented");
    }

    public double getUpperBound(int n) {
        throw new RuntimeException("unimplemented");
    }

    public DemographicFunction getCopy() {
        throw new RuntimeException("unimplemented");
    }
}
