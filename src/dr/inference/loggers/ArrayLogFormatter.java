package dr.inference.loggers;

import java.util.ArrayList;
import java.util.List;

/**
 * @author Alexei Drummond
 */
public class ArrayLogFormatter implements LogFormatter {

    String heading;
    String[] labels = null;
    List<String> lines = new ArrayList<String>();
    List<String[]> valuesList = new ArrayList<String[]>();

    public void startLogging(String title) {
    }

    public void logHeading(String heading) {
        this.heading = heading;
    }

    public void logLine(String line) {
        lines.add(line);
    }

    public void logLabels(String[] labels) {
        if (this.labels == null) {
            this.labels = labels;
        } else throw new RuntimeException("logLabels() method should only be called once!");
    }

    public void logValues(String[] values) {
        valuesList.add(values);
    }

    public void stopLogging() {
    }
}
