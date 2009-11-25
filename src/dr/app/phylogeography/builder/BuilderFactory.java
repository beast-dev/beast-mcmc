package dr.app.phylogeography.builder;

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

    public static int nextCount() {
        return ++count;
    }


    private static int count = 0;
}
