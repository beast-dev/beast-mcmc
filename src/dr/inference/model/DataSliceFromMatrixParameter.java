package dr.inference.model;

import dr.util.Attribute;

public class DataSliceFromMatrixParameter implements Attribute<double[]> {

    private final MatrixParameterInterface parameter;
    private final int slice;
    private final SliceDirection direction;
    private final int dim;

    public DataSliceFromMatrixParameter(MatrixParameterInterface parameter,
                                        int slice,
                                        SliceDirection direction) {
        this.parameter = parameter;
        this.slice = slice;
        this.direction = direction;
        this.dim = direction.getLength(parameter);
    }

    @Override
    public String getAttributeName() {
        return parameter.getParameterName() + "." + (slice + 1) + direction.name();
    }

    @Override
    public double[] getAttributeValue() {
        double[] result = new double[dim];
        for (int i = 0; i < dim; ++i) {
            result[i] = direction.getValue(parameter, slice, i);
        }
        return result;
    }

    public enum SliceDirection {
        ROW_WISE("rows") {
            @Override
            public double getValue(MatrixParameterInterface parameter, int slice, int index) {
                return parameter.getParameterValue(slice, index);
            }

            @Override
            public int getLength(MatrixParameterInterface parameter) {
                return parameter.getColumnDimension();
            }

            @Override
            public int getCount(MatrixParameterInterface parameter) {
                return parameter.getRowDimension();
            }
        },
        COLUMN_WISE("columns") {
            @Override
            public double getValue(MatrixParameterInterface parameter, int slice, int index) {
                return parameter.getParameterValue(index, slice);
            }

            @Override
            public int getLength(MatrixParameterInterface parameter) {
                return parameter.getRowDimension();
            }

            @Override
            public int getCount(MatrixParameterInterface parameter) {
                return parameter.getColumnDimension();
            }
        };

        SliceDirection(String name) {
            this.name = name;
        }

        abstract public double getValue(MatrixParameterInterface parameter, int slice, int index);

        abstract public int getLength(MatrixParameterInterface parameter);

        abstract public int getCount(MatrixParameterInterface parameter);

        public String getName() {
            return name;
        }

        private final String name;
    }
}
