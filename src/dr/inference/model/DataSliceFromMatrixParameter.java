package dr.inference.model;

import dr.util.Attribute;

public class DataSliceFromMatrixParameter extends Parameter.Proxy implements Attribute<double[]>{

    private final MatrixParameterInterface parameter;
    private final int slice;
    private final SliceDirection direction;

    public DataSliceFromMatrixParameter(MatrixParameterInterface parameter,
                                        int slice,
                                        SliceDirection direction) {
        super(makeParameterName(parameter, slice, direction),
                direction.getLength(parameter));
        this.parameter = parameter;
        this.slice = slice;
        this.direction = direction;
    }

    private static String makeParameterName(Parameter parameter, int slice, SliceDirection direction) {
        return parameter.getParameterName() + "." + (slice + 1) + "." + direction.name();
    }

    @Override
    public Bounds<Double> getBounds() {
        return parameter.getBounds();
    }

    @Override
    public double getParameterValue(int dim) {
        return direction.getValue(parameter, slice, dim);
    }

    @Override
    public void setParameterValue(int dim, double value) {
        direction.setValue(parameter, slice, dim, value);
    }

    @Override
    public void setParameterValueQuietly(int dim, double value) {
        direction.setValueQuietly(parameter, slice, dim, value);
    }

    public void fireParameterChangedEvent(int index, Parameter.ChangeType type) {
        super.fireParameterChangedEvent(index, type);
    }

    @Override
    public void setParameterValueNotifyChangedAll(int dim, double value) {
        throw new RuntimeException("Not yet implemented");
    }

    public enum SliceDirection {
        ROW_WISE("rows") {
            @Override
            public double getValue(MatrixParameterInterface parameter, int slice, int index) {
                return parameter.getParameterValue(slice, index);
            }

            @Override
            public void setValue(MatrixParameterInterface parameter, int slice, int index, double value) {
                parameter.setParameterValue(slice, index, value);
            }

            @Override
            public void setValueQuietly(MatrixParameterInterface parameter, int slice, int index, double value) {
                parameter.setParameterValueQuietly(slice, index, value);
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
            public void setValue(MatrixParameterInterface parameter, int slice, int index, double value) {
                parameter.setParameterValue(index, slice, value);
            }

            @Override
            public void setValueQuietly(MatrixParameterInterface parameter, int slice, int index, double value) {
                parameter.setParameterValueQuietly(index, slice, value);
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

        abstract public void setValue(MatrixParameterInterface parameter, int slice, int index, double value);

        abstract public void setValueQuietly(MatrixParameterInterface parameter, int slice, int index, double value);

        abstract public int getLength(MatrixParameterInterface parameter);

        abstract public int getCount(MatrixParameterInterface parameter);

        public String getName() {
            return name;
        }

        private final String name;
    }
}
