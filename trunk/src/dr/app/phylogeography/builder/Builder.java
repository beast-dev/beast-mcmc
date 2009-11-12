package dr.app.phylogeography.builder;

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
    
    String getDescription();

    Layer buildLayer();

    JPanel getEditPanel();    
}
