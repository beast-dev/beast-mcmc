package dr.evomodel.tree;

import dr.evolution.tree.NodeAttributeProvider;
import dr.evolution.tree.Tree;
import dr.evolution.tree.NodeRef;
import dr.evolution.util.Taxon;

/**
 * @author Marc A. Suchard
 * @author Philippe Lemey
 */
public class WanderingTaxonLogger implements NodeAttributeProvider {

    public static final String WANDERER = "wanderingTaxonLogger";
    public static final String RELATIVE = "relative";

    public WanderingTaxonLogger(String name, Taxon taxon, Relative relative) {
        if (name == null) {
            this.name = RELATIVE;
        } else {
            this.name = name;
        }
        this.taxon = taxon;
        this.relative = relative;
    }

    public String[] getNodeAttributeLabel() {
        return new String[] { name };
    }

    public String[] getAttributeForNode(Tree tree, NodeRef node) {
        int rtnValue = 0;
        if (relative == Relative.PARENT) {
            if (isAnyChildEqualToTaxon(tree, node, taxon, null)) {
                rtnValue = 1;
            }
        } else if (relative == Relative.SISTER && !tree.isRoot(node)) {
            if (isAnyChildEqualToTaxon(tree, tree.getParent(node), taxon, node)) {              
                rtnValue = 1;
            }
        }

        return new String[]{Integer.toString(rtnValue)};
    }

    private boolean isAnyChildEqualToTaxon(Tree tree, NodeRef node, Taxon taxon, NodeRef exclude) {
        for (int i = 0; i < tree.getChildCount(node); i++) {
            NodeRef child = tree.getChild(node, i);
            if (child != exclude && tree.isExternal(child)) {
                String taxonString = tree.getNodeTaxon(child).getId();
                if (taxonString.equals(taxon.getId())) {
                    return true;
                }
            }
        }
        return false;
    }

    public enum Relative {
        PARENT,
        SISTER
    }
  
    private String name;
    private Taxon taxon;
    private Relative relative;

}
