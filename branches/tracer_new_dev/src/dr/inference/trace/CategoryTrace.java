package dr.inference.trace;

/**
 *
 */
public class CategoryTrace extends Trace<String> {
        public CategoryTrace(String name, int initialSize) {
            super(name, initialSize);
        }

        public CategoryTrace(String name, String[] values) {
            super(name, values);
        }

        public void add(String value) {
            if (valueCount == values.length) {
                String[] newValues = new String[valueCount + INCREMENT_SIZE];
                System.arraycopy(values, 0, newValues, 0, values.length);
                super.values = newValues;
            }

            super.values[valueCount] = value;
            super.valueCount++;
        }

        public String parserValueWithType(String value) {
            return value;
        }
    }
