package dr.util;

/*
 * ParameterSubsetLogger.java
 *
 * A minimal logger that outputs only a subset of values from a Parameter.
 *
 * Conventions:
 *  - Indices are 0-based internally.
 *  - Column labels are 1-based in the printed name (more human friendly).
 */

import dr.inference.loggers.LogColumn;
import dr.inference.loggers.Loggable;
import dr.inference.loggers.NumberColumn;
import dr.inference.model.Parameter;
import dr.xml.Reportable;

import java.util.Arrays;

public final class ParameterSubsetLogger implements Loggable, Reportable {

    private final Parameter parameter;
    private final int[] indices;          // 0-based indices into the parameter
    private final String prefix;          // column name prefix (defaults to parameter id)

    private LogColumn[] columns;

    /**
     * Log explicit indices (0-based).
     */
    public ParameterSubsetLogger(final Parameter parameter, final int[] indices) {
        this(parameter, indices, null);
    }

    /**
     * Log explicit indices (0-based), with optional name prefix.
     */
    public ParameterSubsetLogger(final Parameter parameter, final int[] indices, final String prefix) {
        if (parameter == null) throw new IllegalArgumentException("parameter cannot be null");
        if (indices == null) throw new IllegalArgumentException("indices cannot be null");
        if (indices.length == 0) throw new IllegalArgumentException("indices cannot be empty");

        this.parameter = parameter;
        this.indices = Arrays.copyOf(indices, indices.length);
        this.prefix = (prefix != null) ? prefix : parameter.getId();

        validateIndices();
    }

    /**
     * Log the first n entries (i = 0,1,...,n-1).
     */
    public ParameterSubsetLogger(final Parameter parameter, final int nFirst) {
        this(parameter, nFirst, null);
    }

    /**
     * Log the first n entries (i = 0,1,...,n-1), with optional name prefix.
     */
    public ParameterSubsetLogger(final Parameter parameter, final int nFirst, final String prefix) {
        if (parameter == null) throw new IllegalArgumentException("parameter cannot be null");
        if (nFirst <= 0) throw new IllegalArgumentException("nFirst must be > 0");

        this.parameter = parameter;
        this.indices = new int[nFirst];
        for (int i = 0; i < nFirst; i++) indices[i] = i;
        this.prefix = (prefix != null) ? prefix : parameter.getId();

        validateIndices();
    }

    private void validateIndices() {
        final int dim = parameter.getDimension();
        for (int idx : indices) {
            if (idx < 0 || idx >= dim) {
                throw new IllegalArgumentException(
                        "Index " + idx + " out of bounds for parameter dimension " + dim);
            }
        }
    }

    @Override
    public LogColumn[] getColumns() {
        if (columns == null) {
            columns = makeColumns();
        }
        return columns;
    }

    private LogColumn[] makeColumns() {
        final LogColumn[] cols = new LogColumn[indices.length];

        for (int k = 0; k < indices.length; k++) {
            final int i = indices[k];
            final String name = prefix + "[" + (i + 1) + "]"; // 1-based label
            cols[k] = new NumberColumn(name) {
                @Override
                public double getDoubleValue() {
                    return parameter.getParameterValue(i);
                }
            };
        }
        return cols;
    }

    @Override
    public String getReport() {
        final StringBuilder sb = new StringBuilder();
        for (LogColumn c : getColumns()) {
            sb.append(c.getFormatted()).append(' ');
        }
        sb.append('\n');
        return sb.toString();
    }
}
