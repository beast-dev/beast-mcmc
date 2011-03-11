package dr.evolution.tree;

import java.util.ArrayList;
import java.util.List;

/**
 * @author Alexei Drummond
 */
public class RankedNode {

    int rank;
    int n;
    RankedNode child1, child2;
    Clade clade;

    public RankedNode(int i, int n) {
        rank = i - n;
        this.n = n;
        clade = new TerminalClade(rank);
        child1 = null;
        child2 = null;
    }

    public RankedNode(int rank, RankedNode child1, RankedNode child2) {

        this.rank = rank;
        n = child1.n;

        if (child1.clade.disjoint(child2.clade)) {

            clade = new Clade(child1.clade, child2.clade);
            this.child1 = child1;
            this.child2 = child2;
        } else throw new IllegalArgumentException();
    }

    class Clade {

        boolean[] in;

        public Clade() {
            in = new boolean[n];
        }

        public Clade(Clade c1, Clade c2) {
            in = new boolean[n];

            for (int i = 0; i < n; i++) {
                if (c1.contains(i - n) && c2.contains(i - n)) {
                    throw new IllegalArgumentException();
                } else {
                    in[i] = c1.contains(i - n) || c2.contains(i - n);
                }
            }
        }

        public void add(int node) {
            in[node + n] = true;
        }

        public boolean contains(int node) {
            return in[node + n];
        }

        public void set(Clade clade) {
            System.arraycopy(clade.in, 0, in, 0, in.length);
        }

        public boolean equals(Clade c) {
            for (int i = 0; i < in.length; i++) {
                if (c.in[i] != in[i]) return false;
            }
            return true;
        }

        public boolean compatible(Clade clade) {

            return outsize(clade) == 0 || insize(clade) == 0 || insize(clade) == clade.size();
        }

        public int size() {
            int size = 0;
            for (int i = 0; i < n; i++) {
                size += in[i] ? 1 : 0;
            }
            return size;
        }

        private int insize(Clade clade) {
            int insize = 0;
            for (int i = 0; i < n; i++) {
                if (contains(i - n)) {
                    if (clade.contains(i - n)) {
                        insize += 1;
                    }
                }
            }
            return insize;
        }

        private int outsize(Clade clade) {
            int outsize = 0;
            for (int i = 0; i < n; i++) {
                if (contains(i - n)) {
                    if (!clade.contains(i - n)) {
                        outsize += 1;
                    }
                }
            }
            return outsize;
        }


        public boolean disjoint(Clade clade) {
            for (int i = 0; i < n; i++) {
                if (contains(i - n) && clade.contains(i - n)) {
                    return false;
                }
            }
            return true;
        }

        public String toString() {
            StringBuilder builder = new StringBuilder();
            for (boolean b : in) {
                if (b) builder.append('1');
                else builder.append('0');
            }
            return builder.toString();
        }

    }

    class TerminalClade extends Clade {
        int terminal;

        public TerminalClade(int terminal) {
            this.terminal = terminal;
        }

        public void add(int node) {
            throw new IllegalArgumentException();
        }

        public boolean contains(int node) {
            return terminal == node;
        }

        public void set(Clade clade) {
            throw new IllegalArgumentException();
        }
    }

    public static void main(String[] args) {

        int n = 6;

        RankedForest start = new RankedForest.Default(n);

        List<RankedNode> complete = new ArrayList<RankedNode>();

        List<Clade> constraints = new ArrayList<Clade>();
        constraints.add(start.getNodes().get(0).createClade("000011"));
        constraints.add(start.getNodes().get(0).createClade("110000"));


        int[] count = new int[]{0};
        processHistory(start, complete, constraints, count);

        System.out.println("n = " + n);

        System.out.println("Constraints:");
        for (Clade constraint : constraints) {
            System.out.println("  " + constraint);
        }
        System.out.println(complete.size() + " histories found in " + count[0] + " calls");
    }

    private Clade createClade(String s) {
        Clade c = new Clade();
        for (int i = 0; i < s.length(); i++) {
            if (s.charAt(i) == '1') {
                c.add(i - n);
            }
        }
        return c;
    }

    private static void processHistory(RankedForest history, List<RankedNode> completeHistories, List<Clade> constraints, int[] count) {

        count[0] += 1;

        List<RankedNode> nodes = history.getNodes();

        if (nodes.size() == 1) {
            completeHistories.add(nodes.get(0));
        } else {

            for (int i = 0; i < nodes.size(); i++) {
                for (int j = i + 1; j < nodes.size(); j++) {

                    RankedNode parent = new RankedNode(history.rank() + 1, nodes.get(i), nodes.get(j));

                    boolean compatible = true;
                    if (constraints != null) {
                        for (Clade constraint : constraints) {
                            if (!parent.clade.compatible(constraint)) {
                                //System.out.println(parent.clade + " not compatible with " + constraint);
                                compatible = false;
                                break;
                            }
                        }
                    }
                    if (compatible) {
                        RankedForest newHis = new RankedForest.Parent(history, parent);

                        // check if order of constraints okay
                        if (newHis.compatible(constraints)) {
                            processHistory(newHis, completeHistories, constraints, count);
                        }
                    }
                }
            }
        }
    }

}
