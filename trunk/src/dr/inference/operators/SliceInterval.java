package dr.inference.operators;

import dr.inference.model.Likelihood;
import dr.inference.model.Variable;
import dr.inference.prior.Prior;
import dr.math.MathUtils;

/**
 * Constructs a univariate slice sampler interval
 *
 * @author Marc Suchard
 */
public interface SliceInterval {

    public double drawFromInterval(Prior prior, Likelihood likelihood, double cutoffDensity, double width);

    public void setSliceSampler(SliceOperator sliceSampler);

    public abstract class Abstract implements SliceInterval {

        public double drawFromInterval(Prior prior, Likelihood likelihood, double cutoffDensity, double width) {
            double x0 = variable.getValue(0);
            Interval interval = constructInterval(prior, likelihood, x0, cutoffDensity, width);
            double x1 = x0;
            boolean found = false;
            while (!found) {
                x1 = MathUtils.uniform(interval.lower, interval.upper);
                if (cutoffDensity < evaluate(prior, likelihood, x1) && test(prior, likelihood, x0, x1,
                        cutoffDensity, width)) {
                    found = true;
                } else {
                    shrinkInterval(interval, x0, x1);
                }
            }
            return x1;
        }

        public abstract Interval constructInterval(Prior prior, Likelihood likelihood, double x0, 
                                                   double cutoffDensity, double width);


        protected abstract boolean test(Prior prior, Likelihood likelihood, double x0, double x1,
                               double cutoffDensity, double width);

        public void shrinkInterval(Interval interval, double x0, double x1) {
            // Taken from Fig 5 in Neal (2003)
            if (x1 < x0) {
                interval.lower = x1;
            } else {
                interval.upper = x1;
            }
        }

        protected class Interval {
            double lower;
            double upper;

            Interval(double lower, double upper) {
                this.lower = lower;
                this.upper = upper;
            }
        }

        public void setSliceSampler(SliceOperator sliceSampler) {
            this.sliceSampler = sliceSampler;
            this.variable = sliceSampler.getVariable();
         }

        protected double evaluate(Prior prior, Likelihood likelihood, double x) {
            variable.setValue(0, x);
            return sliceSampler.evaluate(likelihood, prior, 1.0);
        }

        protected Variable<Double> variable;
        protected SliceOperator sliceSampler;
    }

    public class SteppingOut extends Abstract {

        public SteppingOut() {
            this(10); // TODO Pick better default
        }

        public SteppingOut(int m) {
            this.m = m;
        }

        public Interval constructInterval(Prior prior, Likelihood likelihood, double x0,
                                          double cutoffDensity, double w) {
            // Taken from Fig 3 in Neal (2003)
            double L = x0 -  w * MathUtils.nextDouble();
            double R = L +  w;
            int J = MathUtils.nextInt(m);
            int K = (m - 1) - J;
            while (J > 0 && cutoffDensity < evaluate(prior, likelihood, L) ) {
                L -=  w;
                J--;
            }
            while (K > 0 && cutoffDensity < evaluate(prior, likelihood, R) ) {
                R +=  w;
                K--;
            }
            return new Interval(L,R);
        }

        protected boolean test(Prior prior, Likelihood likelihood, double x0, double x1,
                               double cutoffDensity, double width) {
            return true;
        }

        private int m; // Maximum number of stepping out intervals to explore
    }

    public class Doubling extends Abstract {

        public Interval constructInterval(Prior prior, Likelihood likelihood, double x0,
                                          double cutoffDensity, double width) {
            // TODO
            return new Interval(0,1);
        }

         protected boolean test(Prior prior, Likelihood likelihood, double x0, double x1,
                               double cutoffDensity, double width) {
             // TODO
            return true;
        }
    }
}
