package dr.app.tools.treeannotator;

import dr.math.MathUtils;

import java.util.*;

/**
 * @author Andrew Rambaut
 * @version $
 */
public final class FingerprintCladeKeys {

    private FingerprintCladeKeys() {}

    public static Long getParentKey(Long key1, Long key2) {
        long fingerprint = key1 ^ key2;
        return fingerprint;
    }

    public static Long getTaxonKey(int taxon) {
        long fingerprint = TAXON_FINGERPRINTS.computeIfAbsent(taxon, integer -> {
            long fp = MathUtils.nextLong();
            // abundance of caution - check each fingerprint is unique
//            while (FINGERPRINTS.contains(fp)) {
//                fp = MathUtils.nextLong();
//            }
//            FINGERPRINTS.add(fp);
            return fp;
        });
        TAXON_FINGERPRINTS.put(taxon, fingerprint);
        return fingerprint;
    }

    private final static Map<Integer, Long> TAXON_FINGERPRINTS = new HashMap<>();
    private final static Set<Long> FINGERPRINTS = new HashSet<>();
}
