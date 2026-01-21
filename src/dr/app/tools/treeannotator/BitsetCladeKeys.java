package dr.app.tools.treeannotator;

/**
 * @author Andrew Rambaut
 * @version $
 */
public class BitsetCladeKeys implements CladeKeys {
    public final static CladeKeys INSTANCE = new BitsetCladeKeys();

    @Override
    public Object getParentKey(Object key1, Object key2) {
        int maxIndex;
        if (key1 instanceof Integer) {
            maxIndex = (Integer) key1;
        } else {
            assert key1 instanceof BitsetKey;
            maxIndex = ((BitsetKey) key1).getMaxIndex();
        }
        if (key2 instanceof Integer) {
            maxIndex = Math.max(maxIndex, (Integer) key2);
        } else {
            assert key2 instanceof BitsetKey;
            maxIndex = Math.max(maxIndex, ((BitsetKey) key2).getMaxIndex());
        }

        BitsetKey key = new BitsetKey(maxIndex);
        if (key1 instanceof Integer) {
            key.set((Integer) key1);
        } else {
            key.setTo((BitsetKey) key1);
        }
        if (key2 instanceof Integer) {
            key.set((Integer) key2);
        } else {
            key.or((BitsetKey) key2);
        }

        return key;
    }

    @Override
    public Object getTaxonKey(int taxon) {
        return taxon;
    }
}
