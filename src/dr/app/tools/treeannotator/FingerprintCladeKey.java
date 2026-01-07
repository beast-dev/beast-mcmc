package dr.app.tools.treeannotator;

import dr.math.MathUtils;

import java.util.*;

/**
 * @author Andrew Rambaut
 * @version $
 */
public class FingerprintCladeKey implements CladeKey {
    public FingerprintCladeKey(int taxon) {
        fingerprint = TAXON_FINGERPRINTS.computeIfAbsent(taxon, integer -> {
            long fp = MathUtils.nextLong();
            // abundance of caution - check each fingerprint is unique
//            while (FINGERPRINTS.contains(fp)) {
//                fp = MathUtils.nextLong();
//            }
//            FINGERPRINTS.add(fp);
            return fp;
        });
        KEYS.put(fingerprint, this);
    }

    public FingerprintCladeKey(FingerprintCladeKey key1, FingerprintCladeKey key2) {
        fingerprint = key1.fingerprint ^ key2.fingerprint;
//        FINGERPRINTS.add(fingerprint);
    }

    public static CladeKey getKey(int taxon) {
        return KEYS.get(TAXON_FINGERPRINTS.get(taxon));
    }

    public static CladeKey getParentKey(FingerprintCladeKey key1, FingerprintCladeKey key2) {
        return new FingerprintCladeKey((FingerprintCladeKey)key1, (FingerprintCladeKey)key2);
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (!(o instanceof FingerprintCladeKey)) return false;
        FingerprintCladeKey that = (FingerprintCladeKey) o;
        return fingerprint == that.fingerprint;
    }

    @Override
    public int hashCode() {
        return Objects.hash(fingerprint);
    }

    private final long fingerprint;

    private final static Map<Long, CladeKey> KEYS = new HashMap<>();
    private final static Map<Integer, Long> TAXON_FINGERPRINTS = new HashMap<>();
    private final static Set<Long> FINGERPRINTS = new HashSet<>();
}
