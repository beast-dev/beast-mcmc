package dr.evolution.tree;

import dr.evolution.util.TaxonList;
import dr.evomodel.branchratemodel.LocalClockModel;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.AbstractModel;
import dr.inference.model.Model;
import dr.inference.model.Variable;

import java.util.BitSet;
import java.util.HashSet;
import java.util.Set;

/**
 * @author Marc A. Suchard
 */
public interface TreeNodeFilter {

    boolean includeNode(Tree tree, NodeRef node);

    public class ExternalInternalNodeFilter implements TreeNodeFilter {

        public ExternalInternalNodeFilter(boolean includeExternalNodes, boolean includeInternalNodes) {
            this.includeExternalNodes = includeExternalNodes;
            this.includeInternalNodes = includeInternalNodes;
        }

        public boolean includeNode(Tree tree, NodeRef node) {
            return ((includeExternalNodes && tree.isExternal(node))
                    || (includeInternalNodes && !tree.isExternal(node)));
        }

        private final boolean includeExternalNodes;
        private final boolean includeInternalNodes;
    }
}
