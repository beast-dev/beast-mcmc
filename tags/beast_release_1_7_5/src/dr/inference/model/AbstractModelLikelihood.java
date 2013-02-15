package dr.inference.model;

import dr.inference.loggers.LogColumn;
import dr.inference.loggers.NumberColumn;

/**
 * @author Joseph Heled
 *         Date: 16/04/2009
 */
public abstract class AbstractModelLikelihood extends AbstractModel implements Likelihood {
    /**
     * @param name Model Name
     */
    public AbstractModelLikelihood(String name) {
        super(name);
    }

    public String prettyName() {
        return Likelihood.Abstract.getPrettyName(this);
    }

    @Override
    public boolean isUsed() {
        return isUsed;
    }

    public void setUsed() {
        isUsed = true;
    }

    public boolean evaluateEarly() {
        return false;
    }

    private boolean isUsed = false;

    // **************************************************************
    // Loggable IMPLEMENTATION
    // **************************************************************

    public LogColumn[] getColumns() {
        return new LogColumn[]{
                new LikelihoodColumn(getId())
        };
    }

    protected class LikelihoodColumn extends NumberColumn {
        public LikelihoodColumn(String label) {
            super(label);
        }

        public double getDoubleValue() {
            return getLogLikelihood();
        }
    }
}
