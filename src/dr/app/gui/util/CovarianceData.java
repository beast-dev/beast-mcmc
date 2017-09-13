package dr.app.gui.util;

import java.util.HashMap;
import java.util.List;
import java.util.Set;

/**
 * @author Guy Baele
 */
public class CovarianceData {

    private HashMap<String, List<Double>> data;

    public CovarianceData() {
        data = new HashMap<String, List<Double>>();
    }

    public void add(String traceName, List values) {
        if (data.containsKey(traceName)) {
            List temp = data.get(traceName);
            temp.addAll(values);
            data.replace(traceName, temp);
        } else {
            data.put(traceName, values);
        }
    }

    public Set<String> getTraceNames() {
        return data.keySet();
    }

    public List getDataForKey(String traceName) {
        return data.get(traceName);
    }

    public int numberOfEntries() {
        return data.size();
    }

    public void clear() {
        this.data.clear();
    }

}
