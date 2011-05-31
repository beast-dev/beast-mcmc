package dr.app.beagle.evomodel.utilities;

import dr.app.beagle.evomodel.treelikelihood.MarkovJumpsBeagleTreeLikelihood;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeTrait;
import dr.inference.loggers.LogColumn;
import dr.inference.loggers.Loggable;
import dr.inference.markovjumps.StateHistory;

import java.util.StringTokenizer;

/**
 * A class to conveniently log a complete state history of a continuous-time Markov chain along a tree
 * simulated using the Uniformization Method
 * <p/>
 * This work is supported by NSF grant 0856099
 * <p/>
 * Minin VN and Suchard MA (2008) Counting labeled transitions in continous-time Markov models of evolution.
 * Journal of Mathematical Biology, 56, 391-412.
 * <p/>
 * Rodrigue N, Philippe H and Lartillot N (2006) Uniformization for sampling realizations of Markov processes:
 * applications to Bayesian implementations of codon substitution models. Bioinformatics, 24, 56-62.
 * <p/>
 * Hobolth A and Stone E (2009) Simulation from endpoint-conditioned, continuous-time Markov chains on a finite
 * state space, with applications to molecular evolution. Annals of Applied Statistics, 3, 1204-1231.
 *
 * @author Marc A. Suchard
 * @author Philippe Lemey
 * @author Andrew Rambaut
 */
public class CompleteHistoryLogger implements Loggable {

    public static final String TOTAL_COUNT_NAME = "totalChangeCount";
    public static final String COMPLETE_HISTORY_NAME = "completeHistory";

    public CompleteHistoryLogger(MarkovJumpsBeagleTreeLikelihood treeLikelihood) {
        this.tree = treeLikelihood.getTreeModel();
        treeTraitHistory = treeLikelihood.getTreeTrait(MarkovJumpsBeagleTreeLikelihood.HISTORY);
        treeTraitCount = treeLikelihood.getTreeTrait(MarkovJumpsBeagleTreeLikelihood.TOTAL_COUNTS);
        if (treeTraitHistory == null) {
            throw new RuntimeException("Tree '" + treeLikelihood.getId() + "' does not have a complete history trait.");
        }
        if (treeTraitCount == null) {
            throw new RuntimeException("No sum");
        }
    }

    public LogColumn[] getColumns() {

        LogColumn[] columns = new LogColumn[2];
        columns[0] = new LogColumn.Abstract(TOTAL_COUNT_NAME) {

            @Override
            protected String getFormattedValue() {
                return treeTraitCount.getTraitString(tree, null);
            }
        };
        columns[1] = new LogColumn.Abstract(COMPLETE_HISTORY_NAME) {

            @Override
            protected String getFormattedValue() {
                boolean empty = true;
                StringBuilder bf = new StringBuilder("{");
                int count = 0;
                for (int i = 0; i < tree.getNodeCount(); ++i) {
                    NodeRef node = tree.getNode(i);
                    if (!tree.isRoot(node)) {
                        NodeRef parent = tree.getParent(node);
                        double parentTime = tree.getNodeHeight(parent);
                        double childTime = tree.getNodeHeight(node);
                        double minTime = Math.min(parentTime, childTime);
                        double maxTime = Math.max(parentTime, childTime);
                        String trait = treeTraitHistory.getTraitString(tree,node);
                        if (trait.compareTo("{}") != 0) {
                            trait = trait.substring(1,trait.length()-1);
                            StringTokenizer st = new StringTokenizer(trait,",");
                            while(st.hasMoreTokens()) {
                                if (!empty) {
                                    bf.append(",");
                                }
                                String event = st.nextToken();
                                event = event.substring(1,event.length()-1);
                                StringTokenizer value = new StringTokenizer(event,":");
                                String source = value.nextToken();
                                String dest = value.nextToken();
                                double thisTime = Double.parseDouble(value.nextToken());
                                if (thisTime < 0.0) {
                                    throw new RuntimeException("negative time");
                                }
                                if (thisTime > maxTime || thisTime < minTime) {
                                    throw new RuntimeException("Invalid simulation time");
                                }

                                StateHistory.addEventToStringBuilder(bf, source, dest,
                                        thisTime);
                                count++;
                                empty = false;
                            }
                        }
                    }
                }                
                bf.append("}").append(" ").append(count);
                return bf.toString();
            }
        };
        return columns;
    }

    private Tree tree;
    private TreeTrait treeTraitHistory;
    private TreeTrait treeTraitCount;
}
