package dr.inference.operators;

import dr.math.MathUtils;

/**
 * Constructs a univariate slice sampler interval
 *
 * @author Marc Suchard
 */
public interface SliceInterval {

    public double drawFromInterval(SliceOperator sliceSampler, double cutoffDensity, double width);

    public abstract class Abstract implements SliceInterval {

        public double drawFromInterval(SliceOperator sliceSampler, double sliceDensity, double width) {
            Interval interval = constructInterval(sliceSampler, sliceDensity, width);
            return interval.lower + (interval.upper - interval.lower) * MathUtils.nextDouble();
        }

        public abstract Interval constructInterval(SliceOperator sliceSampler, double sliceDensity, double width);

        protected class Interval {
            double lower;
            double upper;

            Interval(double lower, double upper) {
                this.lower = lower;
                this.upper = upper;
            }
        }
    }

    public class Doubling extends Abstract {

        public Interval constructInterval(SliceOperator sliceSampler, double sliceDensity, double width) {
            // Do nothing
            return new Interval(0,1);
        }
    }

    public class SteppingOut extends Abstract {

        public Interval constructInterval(SliceOperator sliceSampler, double sliceDensity, double width) {
            // Do nothing
            return new Interval(0,1);
        }
    }
}
