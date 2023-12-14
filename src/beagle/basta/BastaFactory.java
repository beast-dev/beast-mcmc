package beagle.basta;

import beagle.*;

import java.util.logging.Logger;

public class BastaFactory extends BeagleFactory {

    public static BeagleBasta loadBastaInstance(
            int tipCount,
            int partialsBufferCount,
            int compactBufferCount,
            int stateCount,
            int patternCount,
            int eigenBufferCount,
            int matrixBufferCount,
            int categoryCount,
            int scaleBufferCount,
            int[] resourceList,
            long preferenceFlags,
            long requirementFlags) {

        getBeagleJNIWrapper();
        if (BeagleJNIWrapper.INSTANCE != null) {

            getBastaJNIWrapper();
            if (BastaJNIWrapper.INSTANCE != null) {

                try {
                    BeagleBasta beagle = new BastaJNIImpl(
                            tipCount,
                            partialsBufferCount,
                            compactBufferCount,
                            stateCount,
                            patternCount,
                            eigenBufferCount,
                            matrixBufferCount,
                            categoryCount,
                            scaleBufferCount,
                            resourceList,
                            preferenceFlags,
                            requirementFlags
                    );

                    // In order to know that it was a CPU instance created, we have to let BEAGLE
                    // to make the instance and then override it...

                    InstanceDetails details = beagle.getDetails();

                    if (details != null) // If resourceList/requirements not met, details == null here
                        return beagle;

                } catch (BeagleException beagleException) {
                    Logger.getLogger("beagle").info("  " + beagleException.getMessage());
                }
            } else {
                throw new RuntimeException("No acceptable BEAGLE-BASTA library plugin found. " +
                        "Make sure that BEAGLE-BASTA is properly installed or try changing resource requirements.");
            }
        }

        throw new RuntimeException("No acceptable BEAGLE library plugins found. " +
                "Make sure that BEAGLE is properly installed or try changing resource requirements.");
    }

    private static BastaJNIWrapper getBastaJNIWrapper() {
        if (BastaJNIWrapper.INSTANCE == null) {
            try {
                BastaJNIWrapper.loadBastaLibrary();
            } catch (UnsatisfiedLinkError ule) {
                System.err.println("Failed to load BEAGLE-BASTA library: " + ule.getMessage());
            }
        }

        return BastaJNIWrapper.INSTANCE;
    }
}
