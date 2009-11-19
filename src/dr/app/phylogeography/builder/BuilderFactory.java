package dr.app.phylogeography.builder;

import dr.app.phylogeography.spread.SpreadDocument;

import javax.swing.*;

/**
 * @author Andrew Rambaut
 * @version $Id$
 */
public abstract class BuilderFactory {
    
    public abstract String getBuilderName();

    public abstract Builder createBuilder();

    public String toString() {
        return getBuilderName();
    }
}
