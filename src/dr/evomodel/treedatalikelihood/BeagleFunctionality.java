package dr.evomodel.treedatalikelihood;

import beagle.BeagleInfo;

import java.util.ArrayList;
import java.util.List;

/**
 * @author Marc A. Suchard
 */
public class BeagleFunctionality {

    private static boolean checkGTEVersion(int[] versionNumbers){
        int[] beagleVersionNumbers = BeagleInfo.getVersionNumbers();
        if (versionNumbers.length == 0 || beagleVersionNumbers.length == 0)
            return false;
        for (int i = 0; i < versionNumbers.length && i < beagleVersionNumbers.length; i++){
            if (beagleVersionNumbers[i] > versionNumbers[i])
                return true;
            if (beagleVersionNumbers[i] < versionNumbers[i])
                return false;
        }
        return true;
    }

    public static boolean IS_THREAD_COUNT_COMPATIBLE() {
        return checkGTEVersion(new int[]{3,1});
    }

    public static boolean IS_ODD_STATE_SSE_FIXED() {
        // SSE for odd state counts fixed in BEAGLE 3.1.3
        return checkGTEVersion(new int[]{3,1,3});
    }

    static boolean IS_PRE_ORDER_SUPPORTED() {
        int[] versionNumbers = BeagleInfo.getVersionNumbers();
        return checkGTEVersion(new int[]{3,2});
    }

    static boolean IS_MULTI_PARTITION_COMPATIBLE() {
        int[] versionNumbers = BeagleInfo.getVersionNumbers();
        return checkGTEVersion(new int[]{3});
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

    public static List<Long> parseSystemPropertyLongArray(String propertyName) {
        List<Long> order = new ArrayList<>();
        String r = System.getProperty(propertyName);
        if (r != null) {
            String[] parts = r.split(",");
            for (String part : parts) {
                try {
                    long n = Long.parseLong(part.trim());
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
