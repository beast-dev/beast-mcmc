package dr.app.phylogeography.generator;

import dr.app.phylogeography.structure.Layer;
import dr.app.phylogeography.structure.TimeLine;

import java.io.IOException;
import java.io.PrintWriter;
import java.util.Collection;

/**
 * @author Andrew Rambaut
 * @version $Id$
 */
public interface Generator {
    public void generate(PrintWriter writer, final TimeLine timeLine, final Collection<Layer> layers) throws IOException;
}
