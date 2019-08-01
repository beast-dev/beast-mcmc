package dr.evomodel.treedatalikelihood;

import beagle.BeagleInfo;

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
}
