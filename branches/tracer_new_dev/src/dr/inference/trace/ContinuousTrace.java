package dr.inference.trace;

/**
 *
 */
public class ContinuousTrace extends Trace<Double> {
    public ContinuousTrace(String name, int initialSize) {
        super(name, initialSize);
    }

    public ContinuousTrace(String name, Double[] values) {
        super(name, values);
    }

    public void add(Double value) {
        if (valueCount == values.length) {
            Double[] newValues = new Double[valueCount + INCREMENT_SIZE];
            System.arraycopy(values, 0, newValues, 0, values.length);
            super.values = newValues;
        }

        super.values[valueCount] = value;
        super.valueCount++;
    }

    public Double parserValueWithType(String value) {
        return Double.parseDouble(value);
    }
}

