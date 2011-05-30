package dr.app.beagle.evomodel.utilities;

import dr.app.beagle.evomodel.treelikelihood.MarkovJumpsBeagleTreeLikelihood;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeTrait;
import dr.inference.loggers.LogColumn;
import dr.inference.loggers.Loggable;
import dr.inference.loggers.NumberColumn;
import dr.inference.markovjumps.StateHistory;

import java.util.StringTokenizer;

/**
 * @author Marc A. Suchard
 * @author Philippe Lemey
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
                        String trait = treeTraitHistory.getTraitString(tree,node);
                        if (trait.compareTo("{}") != 0) {
                            trait = trait.substring(1,trait.length()-1);
//                            System.err.print(trait + " : ");
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
                                double deltaTime = Double.parseDouble(value.nextToken());
                                double thisTime = parentTime - deltaTime;
                                if (thisTime < 0.0) {
                                    throw new RuntimeException("negative time");
                                }
                                StateHistory.addEventToStringBuilder(bf, source, dest, thisTime);
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
