package dr.app.phylogeography.spread;

import dr.app.phylogeography.builder.Builder;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

/**
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @version $Id$
 */
public class SpreadDocument {
    private final List<InputFile> inputFiles = new ArrayList<InputFile>();
    private final List<Builder> layerBuilders = new ArrayList<Builder>();

    public void addTreeFile(InputFile inputFile) {
        inputFiles.add(inputFile);
        fireDataChanged();
    }

    public List<InputFile> getInputFiles() {
        return inputFiles;
    }

    public void addLayerBuilder(Builder builder) {
        layerBuilders.add(builder);
        fireSettingsChanged();
    }

    public List<Builder> getLayerBuilders() {
        return layerBuilders;
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
