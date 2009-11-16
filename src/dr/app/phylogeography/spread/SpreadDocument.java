package dr.app.phylogeography.spread;

import jebl.evolution.trees.Tree;

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
    private final List<Trees> treeFiles = new ArrayList<Trees>();

    public void addTreeFile(Trees trees) {
        treeFiles.add(trees);
        fireDataChanged();
    }

    public List<Trees> getTreeFiles() {
        return treeFiles;
    }

    class Trees {
        Trees(File file, Tree firstTree) {
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

    private void fireDataChanged() {
        for (Listener listener : listeners) {
            listener.dataChanged();
        }
    }

    private void fireSettingsChanged() {
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
