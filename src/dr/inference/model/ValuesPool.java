package dr.inference.model;

/**
 * @author Joseph Heled
 *         Date: 4/09/2009
 */
public class ValuesPool extends DefaultModel {
    private final Variable<Double> pool;
    private final Variable<Integer> selector;

    ValuesPool(Variable pool, Variable selector) {
        this.pool = pool;
        this.selector = selector;

        assert pool.getSize() == selector.getSize();

        pool.addVariableListener(this);
        selector.addVariableListener(this);
    }

    public Double getValue(int index) {
        final int s = selector.getValue(index);
        return pool.getValue(s);
    }

    public Bounds<Double> getBounds() {
      return bounds;
    }

    private final Bounds bounds = new Bounds<Double>() {
        public Double getUpperLimit(int dimension) {
            final int s = selector.getValue(dimension);
            return pool.getBounds().getUpperLimit(s);
        }

        public Double getLowerLimit(int dimension) {
            final int s = selector.getValue(dimension);
            return pool.getBounds().getLowerLimit(s);
        }

        public int getBoundsDimension() {
            return pool.getSize();
        }
    };
}
