package dr.inference.loggers;

import dr.inference.trace.Trace;
import dr.inference.trace.TraceFactory;

import java.util.ArrayList;
import java.util.List;

/**
 * @author Alexei Drummond
 */
public class ArrayLogFormatter implements LogFormatter {

    String heading;
    String[] labels = null;
    List<String> lines = new ArrayList<String>();
    List<Trace> traces = new ArrayList<Trace>();
    boolean echo = false;

    public ArrayLogFormatter(boolean echo) {
        this.echo = echo;
    }

    public void startLogging(String title) {
    }

    public void logHeading(String heading) {
        this.heading = heading;
        echo(heading);
    }

    public void logLine(String line) {
        lines.add(line);
        echo(line);
    }

    public void logLabels(String[] labels) {
        if (this.labels == null) {
            this.labels = labels;
            for (String label : labels) {
                traces.add(new Trace<Double>(label, TraceFactory.TraceType.DOUBLE));
            }
            echo(labels);
        } else throw new RuntimeException("logLabels() method should only be called once!");
    }

    public void logValues(String[] values) {
        for (int i = 0; i < values.length; i++) {
//            Double v = Double.parseDouble(values[i]);
            traces.get(i).add(Double.parseDouble(values[i]));
        }
        echo(values);
    }

    public void stopLogging() {
    }

    public List<Trace> getTraces() {
        return traces;
    }

    private void echo(String s) {
        if (echo) System.out.println(s);
    }

    private void echo(String[] strings) {
        if (echo) {
            for (String s : strings) {
                System.out.print(s + "\t");
            }
            System.out.println();
        }
    }
}
