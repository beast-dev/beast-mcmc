package dr.app.treespace;

import dr.evolution.io.*;
import dr.evolution.tree.Tree;
import dr.evomodel.tree.AbstractCladeImportanceDistribution;

import java.io.FileReader;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

/**
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @version $Id$
 */
public class TreeSpaceDocument {
    private final List<InputFile> inputFiles = new ArrayList<InputFile>();
    private final CladeSystem cladeSystem = new CladeSystem();

    public void addTreeFile(InputFile inputFile) {
        inputFiles.add(inputFile);
        fireDataChanged();
    }

    public List<InputFile> getInputFiles() {
        return inputFiles;
    }

    public CladeSystem getCladeSystem() {
        return cladeSystem;
    }

    public void fireDataChanged() {
        for (Listener listener : listeners) {
            listener.dataChanged();
        }
    }

    public void fireSettingsChanged() {
        for (Listener listener : listeners) {
            listener.settingsChanged();
        }
    }

    public void addListener(Listener listener) {
        listeners.add(listener);
    }

    public interface Listener {
        void dataChanged();
        void settingsChanged();
    }

    private final Set<Listener> listeners = new HashSet<Listener>();
}
