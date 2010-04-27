package dr.inference.trace;

/**
 * @author Alexei Drummond
 * @author Walter Xie
 */
public class TraceFactory {

    public enum TraceType {
        CONTINUOUS("continuous", "C", Double.class),
        DISCRETE("discrete", "D", Integer.class),
        CATEGORY("category", "S", String.class);

        TraceType(String name, String brief, Class type) {
            this.name = name;
            this.brief = brief;
            this.type = type;
        }

        public String toString() {
            return name;
        }

        public String getBrief() {
            return brief;
        }

        public Class getType() {
            return type;
        }

        private final String name;
        private final String brief;
        private final Class type;
    }

    public static Trace createTrace(TraceType traceType, String name, int initialSize) {

//        Double[] d = new Double[10];
//        Double[] t = new Double[10];
//        System.arraycopy(d, 0, t, 0, d.length);

        System.out.println("create trace (" + name + ") with type " + traceType);
        
        switch (traceType) {
            case CONTINUOUS:
                return new Trace<Double>(name, initialSize, Double.valueOf(0));
            case DISCRETE:
                return new Trace<Integer>(name, initialSize, Integer.valueOf(0));
            case CATEGORY:
                return new Trace<String>(name, initialSize, "initial_value");
        }
        throw new IllegalArgumentException("The trace type " + traceType + " is not recognized.");
    }

}

