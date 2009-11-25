package dr.app.phylogeography.builder;

import dr.app.phylogeography.spread.SpreadDocument;
import dr.app.phylogeography.structure.Layer;

import javax.swing.*;

/**
 * @author Andrew Rambaut
 * @version $Id$
 */
public interface Builder {
    String getBuilderName();

    String getName();
    
    void setName(String name);
    
    String getTableCellContent();

    String getToolTipContent();

    String getDescription();

    SpreadDocument.DataFile getDataFile();

    void setDataFile(SpreadDocument.DataFile dataFile);

    JPanel getEditPanel();

    void setFromEditPanel();

    boolean isBuilt();

    void invalidate();

    void build() throws BuildException;

    Layer getLayer();
}
