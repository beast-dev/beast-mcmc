
package dr.evomodel.branchratemodel;

import dr.evolution.tree.NodeRef;

import java.util.Objects;
import java.util.function.Function;

import static org.netlib.lapack.Dlamch.t;

@FunctionalInterface
public interface CartesianNodeMap {
    void apply(int i, NodeRef nodeI, int j, NodeRef node);
}
