package dr.app.phylogeography.builder;

import dr.app.phylogeography.spread.InputFile;

/**
 * @author Andrew Rambaut
 * @version $Id$
 */
public abstract class BuilderFactory {

    public abstract String getBuilderName();

    public abstract InputFile.Type requiresInputType();

    public abstract Builder createBuilder();

    public String toString() {
        return getBuilderName();
    }

    public static int nextCount() {
        return ++count;
    }


    private static int count = 0;
}
