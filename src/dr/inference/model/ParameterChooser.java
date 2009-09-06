package dr.inference.model;

/**
 * @author Joseph Heled
 *         Date: 31/08/2009
 */
public class ParameterChooser extends Variable.Base<Double> implements Variable<Double> {
    private final ValuesPool pool;
    private final int which[];

    ParameterChooser(String name, ValuesPool pool, int which[]) {
        super(name);
        this.pool = pool;
        this.which = which;
    }

    public int getSize() {
        return which.length;
    }

    public Double getValue(int index) {
        return pool.getValue(which[index]);
    }

    public void setValue(int index, Double value) {
        assert false;
    }

    public Double[] getValues() {
        final int size = getSize();
        Double[] copyOfValues = new Double[size];
        for (int i = 0; i < size; i++) {
            copyOfValues[i] = getValue(i) ;
        }
        return copyOfValues;
    }

    public void storeVariableValues() {
        //
    }

    public void restoreVariableValues() {
        //
    }

    public void acceptVariableValues() {
      //
    }

    private final Bounds<Double> bounds = new Bounds<Double>() {
        public Double getUpperLimit(int dimension) {
            return pool.getBounds().getUpperLimit(which[dimension]);
        }

        public Double getLowerLimit(int dimension) {
            return pool.getBounds().getLowerLimit(which[dimension]);
        }

        public int getBoundsDimension() {
            return getSize();
        }
    };

    public Bounds getBounds() {
        return bounds;
    }

    public void addBounds(Bounds bounds) {
       assert false;
    }
}
