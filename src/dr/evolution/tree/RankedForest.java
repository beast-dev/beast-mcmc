package dr.evolution.tree;

import java.util.ArrayList;
import java.util.List;

/**
 * @author Alexei Drummond
 */
public interface RankedForest {

    public List<RankedNode> getNodes();

    public int getSize();

    public int rank();

    public boolean compatible(List<RankedNode.Clade> constraints);

    public class Default implements RankedForest {

        List<RankedNode> nodes;

        public Default(int size) {
            nodes = new ArrayList<RankedNode>(size);
            for (int i = size - 1; i >= 0; i--) {
                nodes.add(new RankedNode(i, size));
            }
        }

        public int getSize() {
            return nodes.get(0).n;
        }

        public int rank() {
            return nodes.get(0).rank;
        }

        public boolean compatible(List<RankedNode.Clade> constraints) {
            return true;
        }

        public List<RankedNode> getNodes() {
            return nodes;
        }
    }

    public class Parent implements RankedForest {

        RankedForest child;
        RankedNode parent;
        List<RankedNode> nodes;

        public Parent(RankedForest child, RankedNode parent) {
            this.child = child;
            this.parent = parent;
            nodes = new ArrayList<RankedNode>();
            nodes.add(parent);
            nodes.addAll(child.getNodes());
            nodes.remove(parent.child1);
            nodes.remove(parent.child2);
        }

        public List<RankedNode> getNodes() {
            return nodes;
        }

        public int getSize() {
            return nodes.get(0).n;
        }

        public int rank() {
            return nodes.get(0).rank;
        }

        public boolean compatible(List<RankedNode.Clade> constraints) {
            int rank = 0;
            for (RankedNode.Clade constraint : constraints) {
                int newRank = getRank(constraint);
                if (newRank < rank) return false;
                rank = newRank;
            }
            return true;
        }

        private int getRank(RankedNode.Clade constraint) {
            //if (constraint.size() > parent.rank) return Integer.MAX_VALUE;

            if (parent.clade.equals(constraint)) {
                return parent.rank;
            } else {
                if (child instanceof Parent) {
                    return ((Parent) child).getRank(constraint);
                }
                return Integer.MAX_VALUE;
            }
        }
    }
}
