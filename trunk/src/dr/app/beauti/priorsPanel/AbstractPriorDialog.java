package dr.app.beauti.priorsPanel;

import dr.app.beauti.options.Parameter;

/**
 * A base class for dialog boxes for setting the prior. Allows code reuse for opening these
 * in PriorsPanel
 * @author Andrew Rambaut
 * @version $Id$
 */
public interface AbstractPriorDialog {
    int showDialog();

    boolean hasInvalidInput(boolean b);

    void getArguments(Parameter parameter);
}
