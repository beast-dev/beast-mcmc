package dr.app.phylogeography.builder;

import dr.app.phylogeography.spread.InputFile;
import dr.app.phylogeography.spread.MultiLineTableCellContent;
import dr.app.phylogeography.structure.Layer;

import javax.swing.*;

/**
 * @author Andrew Rambaut
 * @version $Id$
 */
public interface Builder extends MultiLineTableCellContent {
    String getBuilderName();

    String getName();

    void setName(String name);

    String getDescription();

    InputFile getInputFile();

    void setInputFile(InputFile inputFile);

    JPanel getEditPanel();

    void setFromEditPanel();

    boolean isBuilt();

    void invalidate();

    void build() throws BuildException;

    Layer getLayer();
}
