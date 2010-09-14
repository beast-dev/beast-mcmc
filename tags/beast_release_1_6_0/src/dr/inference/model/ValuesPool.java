package dr.inference.model;

/**
 * @author Joseph Heled
 *         Date: 4/09/2009
 */
public class ValuesPool extends DefaultModel {
    private final Variable<Double> pool;
    private final Variable<Double> selector;
    private final double defaultValue;

    ValuesPool(Variable<Double> pool, Variable<Double> selector, Double defaultValue) {
        this.pool = pool;
        this.selector = selector;
        this.defaultValue = defaultValue;

        assert pool.getSize() == selector.getSize();

        //pool.addVariableListener(this);
       // selector.addVariableListener(this);

        addVariable(pool);
        addVariable(selector);
        
        addStatistic(numberOfParams);
    }

    public Variable<Double> getPool() {
        return pool;
    }

    public Variable<Double> getSelector() {
        return selector;
    }

    public int length() {
        return pool.getSize();
    }

    private int get(int index) {
       final double s = selector.getValue(index);
       return s < 0 ? -1 : (int)(s+0.5);
    }

    public Double getValue(int index) {
       final int s = get(index);
       return s < 0 ? defaultValue : pool.getValue(s);
    }

    public Bounds<Double> getBounds() {
      return bounds;
    }

    private final Bounds<Double> bounds = new Bounds<Double>() {
        public Double getUpperLimit(int dimension) {
            final int s = get(dimension);
            return s < 0 ? defaultValue : pool.getBounds().getUpperLimit(s);
        }

        public Double getLowerLimit(int dimension) {
            final int s = get(dimension);
            return s < 0 ? defaultValue : pool.getBounds().getLowerLimit(s);
        }

        public int getBoundsDimension() {
            return pool.getSize();
        }
    };

    @SuppressWarnings({"FieldCanBeLocal"})
    private final Statistic numberOfParams = new Statistic.Abstract() {

        public String getStatisticName() {
            return "numberOfParams";
        }

        public int getDimension() {
            return 1;
        }

        public double getStatisticValue(int dim) {
           double m = -1; 
           for(int k = 0; k < selector.getSize(); ++k) {
               final Double sk = selector.getValue(k);
               if( sk > m ) {
                  m = sk;
              }
           }
            return m+1;
        }
    };

    public boolean hasChanged(int i, Object object, int index) {
        if( object == pool ) {
            if( get(i) == index ) {
                return true;
            }
        } else if ( object == selector ) {
            if( i == index ) {
                return true;
            }
        }
        return false;
    }
}
