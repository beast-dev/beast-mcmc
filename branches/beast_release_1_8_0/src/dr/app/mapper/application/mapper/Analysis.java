package dr.app.mapper.application.mapper;

import dr.app.mapper.application.MapperDocument;

import java.util.ArrayList;
import java.util.List;

/**
 * Date: 09/03/2013
 * Time: 09:16
 *
 * @author rambaut
 */
public class Analysis {
    public void addDataSet(MapperDocument dataSet) {
        dataSets.add(dataSet);
        fireDataSetChanged();
    }

    public int getDataSetCount() {
        return dataSets.size();
    }

    public MapperDocument getDataSet(int index) {
        return dataSets.get(index);
    }

    /**
     * Return a copy of the list of data sets.
     * @return the data sets
     */
    public List<MapperDocument> getDataSets() {
        return new ArrayList<MapperDocument>(dataSets);
    }

    // Listeners and broadcasting
    public void addListener(AnalysisListener listener) {
        listeners.add(listener);
    }

    public void removeListener(AnalysisListener listener) {
        listeners.remove(listener);
    }

    private void fireDataSetChanged() {
      for (AnalysisListener listener : listeners) {
          listener.analysisChanged();
      }
    }

    private final List<AnalysisListener> listeners = new ArrayList<AnalysisListener>();
    private final List<MapperDocument> dataSets = new ArrayList<MapperDocument>();
}
