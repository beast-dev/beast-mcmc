package dr.app.beauti.options;

import java.util.List;

/**
 * @author Andrew Rambaut
 * @version $Id$
 */
public interface ComponentOptions {

    void createParameters(ModelOptions modelOptions);

    void selectParameters(ModelOptions modelOptions, List<Parameter> params);

    void selectStatistics(ModelOptions modelOptions, List<Parameter> stats);

    void selectOperators(ModelOptions modelOptions, List<Operator> ops);
}
