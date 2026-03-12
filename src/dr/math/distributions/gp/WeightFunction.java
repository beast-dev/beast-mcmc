package dr.math.distributions.gp;

import java.util.HashMap;
import java.util.Map;
import java.util.function.Supplier;

public interface WeightFunction {
    double getWeight(double x);

    default void configure(Map<String, Double> params) {};

    WeightFunction IDENTITY = x -> 1.0;

    class SigmoidWeightFunction implements WeightFunction {
        private double location = 0.0; // default
        private double scale = 1.0; // default

        @Override
        public void configure(Map<String, Double> params) {
            if (params.containsKey("scale")) {
                scale = params.get("scale");
            }
            if (params.containsKey("location")) {
                location = params.get("location");
            }
        }

        @Override
        public double getWeight(double x) {
            return 1.0 / (1.0 + Math.exp(- scale * ( x - location)));
        }
    }

    class SigmoidComplementWeightFunction implements WeightFunction {
        private double location = 0.0; // default
        private double scale = 1.0; // default

        @Override
        public void configure(Map<String, Double> params) {
            if (params.containsKey("scale")) {
                scale = params.get("scale");
            }
            if (params.containsKey("location")) {
                location = params.get("location");
            }
        }

        @Override
        public double getWeight(double x) {
            return 1.0 - 1.0 / (1.0 + Math.exp(- scale * (x - location)));
        }
    }

    class LinearWeightFunction implements WeightFunction {
        private double slope = 1.0;
        private double intercept = 0.0;

        @Override
        public void configure(Map<String, Double> params) {
            if (params.containsKey("slope")) slope = params.get("slope");
            if (params.containsKey("intercept")) intercept = params.get("intercept");
        }

        @Override
        public double getWeight(double x) {
            return slope * x + intercept;
        }
    }

    class WeightFunctionFactory {

        private static final Map<String, Supplier<WeightFunction>> registry = new HashMap<>();

        static {
            register("sigmoid", SigmoidWeightFunction::new);
            register("sigmoidComplement", SigmoidComplementWeightFunction::new);
            register("identity", () -> WeightFunction.IDENTITY);
            register("linear", LinearWeightFunction::new);
        }

        public static void register(String name, Supplier<WeightFunction> constructor) {
            registry.put(name.toLowerCase(), constructor);
        }

        public static WeightFunction create(String type, Map<String, Double> params) {
            Supplier<WeightFunction> constructor = registry.get(type.toLowerCase());
            if (constructor == null) {
                throw new IllegalArgumentException("Unknown weight function: " + type);
            }
            WeightFunction wf = constructor.get();
            wf.configure(params);
            return wf;
        }
    }
}