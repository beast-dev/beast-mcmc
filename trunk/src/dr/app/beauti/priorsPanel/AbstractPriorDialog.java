package dr.app.beauti.priorsPanel;

import dr.app.beauti.options.Parameter;

/**
 * @author Andrew Rambaut
 * @version $Id$
 */
public interface AbstractPriorDialog {
    int showDialog();

    boolean hasInvalidInput(boolean b);

    void getArguments(Parameter parameter);
}
