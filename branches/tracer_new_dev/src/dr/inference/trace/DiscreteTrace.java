package dr.inference.trace;

/**
 *
 */
public class DiscreteTrace extends Trace<Integer> {
    public DiscreteTrace(String name, int initialSize) {
        super(name, initialSize);
    }

    public DiscreteTrace(String name, Integer[] values) {
        super(name, values);
    }

    public void add(Integer value) {
        if (valueCount == values.length) {
            Integer[] newValues = new Integer[valueCount + INCREMENT_SIZE];
            System.arraycopy(values, 0, newValues, 0, values.length);
            super.values = newValues;
        }

        super.values[valueCount] = value;
        super.valueCount++;
    }

    public Integer parserValueWithType(String value) {
        return Integer.parseInt(value);
    }
}

