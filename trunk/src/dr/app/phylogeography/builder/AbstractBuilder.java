package dr.app.phylogeography.builder;

import dr.app.phylogeography.structure.*;
import dr.app.phylogeography.spread.SpreadDocument;
import jebl.evolution.trees.RootedTree;
import jebl.evolution.graphs.Node;

import javax.swing.*;
import java.awt.*;

import jam.panels.OptionsPanel;
import org.virion.jam.components.WholeNumberField;
import org.virion.jam.components.RealNumberField;

/**
 * @author Andrew Rambaut
 * @version $Id$
 */
public abstract class AbstractBuilder implements Builder {

    private SpreadDocument.DataFile dataFile;

    public String getName() {
        return name;
    }

    public void setName(final String name) {
        this.name = name;
    }

    public String getDescription() {
        return description;
    }

    public void setDescription(final String description) {
        this.description = description;
    }

    public SpreadDocument.DataFile getDataFile() {
        return dataFile;
    }

    public void setDataFile(SpreadDocument.DataFile dataFile) {
        this.dataFile = dataFile;
    }

    public boolean isVisible() {
        return isVisible;
    }

    public void setVisible(final boolean visible) {
        isVisible = visible;
    }


    private String name;
    private String description;
    private boolean isVisible;
}