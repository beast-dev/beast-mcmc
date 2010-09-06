package dr.app.phylogeography.spread;

import javax.swing.*;
import java.awt.*;

/**
 * @author Andrew Rambaut
 * @version $Id$
 */
public interface DeleteActionResponder {
    Action getDeleteAction();
    void delete();
}
