package dr.app.phylogeography.spread;

import javax.swing.*;

/**
 * @author Andrew Rambaut
 * @version $Id$
 */
public interface MultiLineTableCellContent {

    Icon getTableCellIcon();

    String getTableCellContent();

    String getToolTipContent();
}
