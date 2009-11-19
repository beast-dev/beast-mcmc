package dr.app.phylogeography.builder;

import dr.app.phylogeography.structure.Layer;
import dr.app.phylogeography.spread.SpreadDocument;

import javax.swing.*;

/**
 * @author Andrew Rambaut
 * @version $Id$
 */
public interface Builder {
    String getBuilderName();

    String getName();
    
    void setName(String name);
    
    String getDescription();

    SpreadDocument.DataFile getDataFile();

    void setDataFile(SpreadDocument.DataFile dataFile);

    Layer buildLayer();

    JPanel getEditPanel();

    void setFromEditPanel();
}
