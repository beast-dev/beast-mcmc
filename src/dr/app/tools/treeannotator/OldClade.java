package dr.app.tools.treeannotator;

import dr.evolution.util.Taxon;
import dr.util.Pair;

import java.util.BitSet;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

/**
 * @author Andrew Rambaut
 * @version $
 */
class OldClade implements Clade {
    public OldClade(BitSet bits) {
        this.bits = bits;
        count = 0;
        credibility = 0.0;
        size = bits.cardinality();
    }

    public int getCount() {
        return count;
    }

    public void setCount(int count) {
        this.count = count;
    }

    public double getCredibility() {
        return credibility;
    }

    public void setCredibility(double credibility) {
        this.credibility = credibility;
    }

    @Override
    public int getSize() {
        return size;
    }

    @Override
    public Taxon getTaxon() {
        return taxon;
    }

    @Override
    public Clade getBestLeft() {
        return bestLeft;
    }

    @Override
    public Clade getBestRight() {
        return bestRight;
    }

    public Integer getHash() {
        return hashCode();
    }

    public Integer getKey() {
        return hashCode();
    }

    @Override
    public void addAttributeValues(Object[] values) {

    }

    @Override
    public List<Object[]> getAttributeValues() {
        return attributeValues;
    }

    public void addSubclades(BitSet subClade1, BitSet subClade2) {
        if (this.subClades == null) {
            this.subClades = new HashSet<>();
        }
        // Store the subclade with lowest first set bit index as the first of the pair to make
        // sure the order is the same if the pair is the same.
        if (subClade1.nextSetBit(0) < subClade2.nextSetBit(0)) {
            this.subClades.add(new Pair<>(subClade1, subClade2));
        } else {
            this.subClades.add(new Pair<>(subClade2, subClade1));
        }
    }

    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        final OldClade clade = (OldClade) o;

        return !(bits != null ? !bits.equals(clade.bits) : clade.bits != null);

    }

    public int hashCode() {
        return (bits != null ? bits.hashCode() : 0);
    }

    public String toString() {
        return "clade " + bits.toString();
    }

    int count;
    double credibility;
    final int size;
    final BitSet bits;
    Taxon taxon = null;
    List<Object[]> attributeValues = null;
    Set<Pair<BitSet, BitSet>> subClades = null;
    OldClade bestLeft = null;
    OldClade bestRight = null;
    double bestSubTreeCredibility;
}
