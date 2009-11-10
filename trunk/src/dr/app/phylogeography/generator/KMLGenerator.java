package dr.app.phylogeography.generator;

import dr.app.phylogeography.structure.Layer;

import java.io.PrintWriter;
import java.io.IOException;
import java.util.Collection;

/**
 * @author Andrew Rambaut
 * @version $Id$
 */
public class KMLGenerator implements Generator {
    public void generate(PrintWriter writer, Collection<Layer> layers) throws IOException {
    }

    public String toString() {
        return "KML";
    }
}
