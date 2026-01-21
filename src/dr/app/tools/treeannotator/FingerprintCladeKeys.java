package dr.app.tools.treeannotator;

import dr.math.MathUtils;

import java.util.*;

/**
 * Implements clade keys using the random tip numbers and XOR compositing method of Patrick Varilly
 * described in https://doi.org/10.1101/2025.03.25.645253
 *
 * @author Andrew Rambaut
 * @version $
 */
public final class FingerprintCladeKeys implements CladeKeys {
    public final static CladeKeys INSTANCE = new FingerprintCladeKeys();

    private FingerprintCladeKeys() {}

    @Override
    public Object getParentKey(Object key1, Object key2) {
        long fingerprint = ((Long)key1) ^ ((Long)key2);
        return fingerprint;
    }

    @Override
    public Object getTaxonKey(int taxon) {
        long fingerprint = TAXON_FINGERPRINTS.computeIfAbsent(taxon, integer -> {
            long fp = MathUtils.nextLong();
            // abundance of caution - check each fingerprint is unique
            while (FINGERPRINTS.contains(fp)) {
                assert true : "Taxon fingerprint collision"; // if exceptions are on then flag and stop
                fp = MathUtils.nextLong();
            }
            FINGERPRINTS.add(fp);
            return fp;
        });
        TAXON_FINGERPRINTS.put(taxon, fingerprint);
        return fingerprint;
    }

    private final static Map<Integer, Long> TAXON_FINGERPRINTS = new HashMap<>();
    private final static Set<Long> FINGERPRINTS = new HashSet<>();
}
