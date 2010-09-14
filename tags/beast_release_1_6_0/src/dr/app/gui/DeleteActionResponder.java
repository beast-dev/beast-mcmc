package dr.app.gui;

import javax.swing.*;

/**
 * @author Andrew Rambaut
 * @version $Id$
 */
public interface DeleteActionResponder {
    Action getDeleteAction();
    void delete();
}
