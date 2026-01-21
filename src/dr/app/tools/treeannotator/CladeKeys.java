package dr.app.tools.treeannotator;

/**
 * @author Andrew Rambaut
 * @version $
 */
public interface CladeKeys {
    Object getParentKey(Object key1, Object key2);

    Object getTaxonKey(int taxon);
}
