
package dr.evomodel.branchratemodel;

import dr.evolution.tree.NodeRef;

@FunctionalInterface
public interface NodeRateMap {
    double apply(int index, NodeRef node, double rate);
}
