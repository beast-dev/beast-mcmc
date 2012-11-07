package dr.app.beagle.evomodel.branchmodel;

import dr.evolution.tree.NodeRef;

/**
 * @author Andrew Rambaut
 * @author Filip Bielejec
 * @version $Id$
 */
public class HomogeneousBranchModel implements BranchSpecificModel {
    public Mapping getBranchModelMapping(NodeRef node) {
        return DEFAULT;
    }
}
