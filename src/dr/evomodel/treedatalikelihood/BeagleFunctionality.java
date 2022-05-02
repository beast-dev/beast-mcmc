package dr.evomodel.treedatalikelihood;

import beagle.BeagleInfo;

import java.util.ArrayList;
import java.util.List;

/**
 * @author Marc A. Suchard
 */
public class BeagleFunctionality {

    public static boolean IS_THREAD_COUNT_COMPATIBLE() {
        int[] versionNumbers = BeagleInfo.getVersionNumbers();
        return versionNumbers.length != 0 && versionNumbers[0] >= 3 && versionNumbers[1] >= 1;
    }

    public static boolean IS_ODD_STATE_SSE_FIXED() {
        // SSE for odd state counts fixed in BEAGLE 3.1.3
        int[] versionNumbers = BeagleInfo.getVersionNumbers();
        return versionNumbers.length != 0 && versionNumbers[0] >= 3 && versionNumbers[1] >= 1 && versionNumbers[2] >= 3;
    }

    static boolean IS_PRE_ORDER_SUPPORTED() {
        int[] versionNumbers = BeagleInfo.getVersionNumbers();
        return versionNumbers.length != 0 && versionNumbers[0] >= 3 && versionNumbers[1] >= 2;
    }

    static boolean IS_MULTI_PARTITION_COMPATIBLE() {
        int[] versionNumbers = BeagleInfo.getVersionNumbers();
        return versionNumbers.length != 0 && versionNumbers[0] >= 3;
    }

    public static List<Integer> parseSystemPropertyIntegerArray(String propertyName) {
        List<Integer> order = new ArrayList<>();
        String r = System.getProperty(propertyName);
        if (r != null) {
            String[] parts = r.split(",");
            for (String part : parts) {
                try {
                    int n = Integer.parseInt(part.trim());
                    order.add(n);
                } catch (NumberFormatException nfe) {
                    System.err.println("Invalid entry '" + part + "' in " + propertyName);
                }
            }
        }
        return order;
    }

    public static List<String> parseSystemPropertyStringArray(String propertyName) {

        List<String> order = new ArrayList<>();

        String r = System.getProperty(propertyName);
        if (r != null) {
            String[] parts = r.split(",");
            for (String part : parts) {
                try {
                    String s = part.trim();
                    order.add(s);
                } catch (NumberFormatException nfe) {
                    System.err.println("Invalid entry '" + part + "' in " + propertyName);
                }
            }
        }
        return order;
    }
}
