package dr.inference.model;

import dr.inference.loggers.LogColumn;
import dr.inference.loggers.NumberColumn;

/**
 * @author joseph
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

    // **************************************************************
    // Loggable IMPLEMENTATION
    // **************************************************************

    public LogColumn[] getColumns() {
        return new LogColumn[]{
                new NumberColumn(getId()) {
                    public double getDoubleValue() {
                        return getLogLikelihood();
                    }
                }
        };
    }
}
