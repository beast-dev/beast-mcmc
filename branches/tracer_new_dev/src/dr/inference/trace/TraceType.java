package dr.inference.trace;

/**
 * @author Alexei Drummond
 * @author Walter Xie
 */
public interface TraceType<T> {

    public static final int INITIAL_SIZE = 1000;
    public static final int INCREMENT_SIZE = 1000;

    public T parserValueWithType(String value);

    /**
     * @param value the valued to be added
     */
    public void add(T value);

    enum Type {
        CONTINUOUS("continuous", Double.class),
        DISCRETE("discrete", Integer.class),
        CATEGORY("category", String.class),;

        Type(String name, Class type) {
            this.name = name;
            this.type = type;
        }

        public String toString() {
            return name;
        }

        public Class getType() {
            return type;
        }

        private final String name;
        private final Class type;
    }
}

