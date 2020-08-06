package dr.evomodel.tree;

import dr.evolution.tree.*;
import dr.evolution.util.MutableTaxonListListener;
import dr.evolution.util.Taxon;
import dr.inference.model.AbstractModel;
import dr.util.*;
import org.w3c.dom.Document;
import org.w3c.dom.Element;

import java.util.*;

/**
 * @author Andrew Rambaut
 * @version $
 */
public abstract class AbstractTreeModel extends AbstractModel implements MutableTreeModel, Keywordable, Citable {
    private final boolean isVariable;
    private final List<String> keywords = new ArrayList<String>();
    private String id = null;
    private AttributeHelper treeAttributes = null;
    /**
     * holds the units of the trees branches.
     */
    private Type units = Type.SUBSTITUTIONS;

    public AbstractTreeModel(String name, boolean isVariable) {
        super(name);
        this.isVariable = isVariable;
    }

    /**
     * Return the units that this tree is expressed in.
     */
    public Type getUnits() {
        return units;
    }

    /**
     * Sets the units that this tree is expressed in.
     */
    public void setUnits(Type units) {
        this.units = units;
    }

    /**
     * @return a count of the number of taxa in the list.
     */
    public int getTaxonCount() {
        return getExternalNodeCount();
    }

    /**
     * @return the ID of the taxon of the ith external node. If it doesn't have
     *         a taxon, returns the ID of the node itself.
     */
    public String getTaxonId(int taxonIndex) {
        Taxon taxon = getTaxon(taxonIndex);
        if (taxon != null) {
            return taxon.getId();
        } else {
            return null;
        }
    }

    /**
     * returns the index of the taxon with the given id.
     */
    public int getTaxonIndex(String id) {
        for (int i = 0, n = getTaxonCount(); i < n; i++) {
            if (getTaxonId(i).equals(id)) return i;
        }
        return -1;
    }

    /**
     * returns the index of the given taxon.
     */
    public int getTaxonIndex(Taxon taxon) {
        for (int i = 0, n = getTaxonCount(); i < n; i++) {
            if (getTaxon(i) == taxon) return i;
        }
        return -1;
    }

    public List<Taxon> asList() {
        List<Taxon> taxa = new ArrayList<Taxon>();
        for (int i = 0, n = getTaxonCount(); i < n; i++) {
            taxa.add(getTaxon(i));
        }
        return taxa;
    }

    public Iterator<Taxon> iterator() {
        return new Iterator<Taxon>() {
            private int index = -1;

            public boolean hasNext() {
                return index < getTaxonCount() - 1;
            }

            public Taxon next() {
                index++;
                return getTaxon(index);
            }

            public void remove() { /* do nothing */ }
        };
    }

    /**
     * @param taxonIndex the index of the taxon whose attribute is being fetched.
     * @param name       the name of the attribute of interest.
     * @return an object representing the named attributed for the taxon of the given
     *         external node. If the node doesn't have a taxon then the nodes own attribute
     *         is returned.
     */
    public Object getTaxonAttribute(int taxonIndex, String name) {
        Taxon taxon = getTaxon(taxonIndex);
        if (taxon != null) {
            return taxon.getAttribute(name);
        }
        return null;
    }

    public int addTaxon(Taxon taxon) {
        throw new IllegalArgumentException("Cannot add taxon to a TreeModel");
    }

    public boolean removeTaxon(Taxon taxon) {
        throw new IllegalArgumentException("Cannot add taxon to a TreeModel");
    }

    public void setTaxonId(int taxonIndex, String id) {
        throw new IllegalArgumentException("Cannot set taxon id in a TreeModel");
    }

    public void setTaxonAttribute(int taxonIndex, String name, Object value) {
        throw new IllegalArgumentException("Cannot set taxon attribute in a TreeModel");
    }

    public void addMutableTreeListener(MutableTreeListener listener) {
    } // Do nothing at the moment

    public void addMutableTaxonListListener(MutableTaxonListListener listener) {
    } // Do nothing at the moment

    /**
     * @return the id.
     */
    public String getId() {
        return id;
    }

    /**
     * Sets the id.
     */
    public void setId(String id) {
        this.id = id;
    }

    public Element createElement(Document document) {
        throw new RuntimeException("Not implemented yet");
    }

    /**
     * Sets an named attribute for this object.
     *
     * @param name  the name of the attribute.
     * @param value the new value of the attribute.
     */
    public void setAttribute(String name, Object value) {
        if (treeAttributes == null)
            treeAttributes = new AttributeHelper();
        treeAttributes.setAttribute(name, value);
    }

    /**
     * @param name the name of the attribute of interest.
     * @return an object representing the named attributed for this object.
     */
    public Object getAttribute(String name) {
        if (treeAttributes == null)
            return null;
        else
            return treeAttributes.getAttribute(name);
    }

    /**
     * @return an iterator of the attributes that this object has.
     */
    public Iterator<String> getAttributeNames() {
        if (treeAttributes == null)
            return null;
        else
            return treeAttributes.getAttributeNames();
    }

    /**
     * @return a string containing a newick representation of the tree
     */
    public final String getNewick() {
        return TreeUtils.newick(this);
    }

    /**
     * @return a string containing a newick representation of the tree
     */
    public String toString() {
        return getNewick();
    }

    public Tree getCopy() {
        throw new UnsupportedOperationException("please don't call this function");
    }

    @Override
    public void addKeyword(String keyword) {
        keywords.add(keyword);
    }

    @Override
    public List<String> getKeywords() {
        return keywords;
    }

    public boolean hasNodeHeights() {
        return true;
    }

    public boolean hasBranchLengths() {
        return true;
    }

    public double getBranchLength(NodeRef node) {
        NodeRef parent = getParent(node);
        if (parent == null) {
            return 0.0;
        }

        return getNodeHeight(parent) - getNodeHeight(node);
    }

    public boolean isTipDateSampled() {
        return false;
    }

    @Override
    public boolean isVariable() {
        return isVariable;
    }

    @Override
    public Citation.Category getCategory() {
        return Citation.Category.TREE_PRIORS;
    }

    @Override
    public String getDescription() {
        return "Sampling tip dates model";
    }

    @Override
    public List<Citation> getCitations() {
        if (isTipDateSampled()) {
            return Arrays.asList(new Citation(
                            new Author[]{
                                    new Author("B", "Shapiro"),
                                    new Author("SYW", "Ho"),
                                    new Author("AJ", "Drummond"),
                                    new Author("MA", "Suchard"),
                                    new Author("OG", "Pybus"),
                                    new Author("A", "Rambaut"),
                            },
                            "A Bayesian phylogenetic method to estimate unknown sequence ages",
                            2010,
                            "Mol Biol Evol",
                            28,
                            879, 887,
                            "10.1093/molbev/msq262"
                    ),
                    new Citation(
                            new Author[]{
                                    new Author("AJ", "Drummond"),
                            },
                            "PhD Thesis",
                            2002,
                            "University of Auckland",
                            ""
                    ));
        } else {
            return Collections.emptyList();
        }
    }

}
