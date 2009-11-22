package dr.app.phylogeography.builder;

import dr.app.phylogeography.spread.SpreadDocument;
import dr.app.phylogeography.structure.Layer;

/**
 * @author Andrew Rambaut
 * @version $Id$
 */
public abstract class AbstractBuilder implements Builder {

    private SpreadDocument.DataFile dataFile = null;

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

    public boolean isBuilt() {
        return layer != null;
    }

    public void invalidate() {
        layer = null;
    }

    public void build() throws BuildException {
        layer = buildLayer();
    }

    protected abstract Layer buildLayer() throws BuildException;

    public Layer getLayer() {
        if (layer == null) {
            throw new IllegalArgumentException("getLayer accessed but layer has not been built");
        }
        return layer;
    }

    private String name;
    private String description;
    private boolean isVisible;

    private Layer layer = null;
}