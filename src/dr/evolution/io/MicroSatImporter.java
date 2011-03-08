package dr.evolution.io;

import java.io.Reader;
import java.io.Writer;

/**
 * @author Alexei Drummond
 * @author Walter Xie
 */
public class MicroSatImporter extends Importer {

    public MicroSatImporter(Reader reader) {
        this(reader, null);
    }

    public MicroSatImporter(Reader reader, Writer commentWriter) {
        super(reader, commentWriter);
    }






}
