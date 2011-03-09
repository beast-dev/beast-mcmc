package dr.evolution.io;

import dr.evolution.alignment.Patterns;

import java.io.IOException;
import java.util.List;

/**
 * @author Alexei Drummond
 * @author Walter Xie
 */
public interface PatternImporter {

    // importPatternList
//	PatternList importPatternList() throws IOException, Importer.ImportException;

    // importPatterns
	List<Patterns> importPatterns() throws IOException, Importer.ImportException;
}
