package dr.app.phylogeography.spread;

import dr.evolution.tree.Tree;
import dr.app.phylogeography.builder.Builder;

import java.util.Set;
import java.util.HashSet;
import java.util.List;
import java.util.ArrayList;
import java.io.File;

/**
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @version $Id$
 */
public class SpreadDocument {
    private final List<DataFile> dataFiles = new ArrayList<DataFile>();
    private final List<Builder> layerBuilders = new ArrayList<Builder>();


    public void addTreeFile(DataFile dataFile) {
        dataFiles.add(dataFile);
        fireDataChanged();
    }

    public List<DataFile> getDataFiles() {
        return dataFiles;
    }

    public void addLayerBuilder(Builder builder) {
        layerBuilders.add(builder);
        fireSettingsChanged();
    }

    public List<Builder> getLayerBuilders() {
        return layerBuilders;
    }

    public static class DataFile {
        DataFile(File file, Tree firstTree) {
            this.file = file;
            this.firstTree = firstTree;
        }

        public File getFile() {
            return file;
        }

        public Tree getFirstTree() {
            return firstTree;
        }

        private final File file;
        private final Tree firstTree;                                
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

    interface Listener {
        void dataChanged();
        void settingsChanged();
    }

    private final Set<Listener> listeners = new HashSet<Listener>();
}
