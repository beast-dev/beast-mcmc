package dr.util;

import java.util.HashMap;
import java.util.Map;

/**
 * @author Marc A. Suchard
 */
public class BenchmarkTimer {

    public void startTimer(String key) {
        startTimes.put(key, getTime());
    }

    public void stopTimer(String key) {

        long end = getTime();
        long start = startTimes.get(key);

        Long total = totals.get(key);
        if (total == null) {
            total = 0L;
        }

        totals.put(key, total + (end - start));
    }
    
    public String toString() {

        StringBuilder sb = new StringBuilder();

        sb.append("\nTIMING:");
        for (String key : totals.keySet()) {
            String value = String.format("%4.3e", (double) totals.get(key));
            sb.append("\n").append(key).append(delimiter).append(value);
        }
        sb.append("\n");

        return sb.toString();
    }

    private long getTime() {
//        return System.nanoTime();
        return System.currentTimeMillis();
    }

    private final Map<String, Long> startTimes = new HashMap<>();
    private final Map<String, Long> totals = new HashMap<>();

    private final static String delimiter = "\t\t";
}
