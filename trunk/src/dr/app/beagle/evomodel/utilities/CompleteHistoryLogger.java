package dr.app.beagle.evomodel.utilities;

import dr.app.beagle.evomodel.treelikelihood.MarkovJumpsBeagleTreeLikelihood;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeTrait;
import dr.inference.loggers.LogColumn;
import dr.inference.loggers.Loggable;
import dr.inference.markovjumps.StateHistory;
import dr.util.Citable;
import dr.util.Citation;
import dr.util.CommonCitations;

import java.util.ArrayList;
import java.util.List;
import java.util.StringTokenizer;
import java.util.logging.Logger;

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
public class CompleteHistoryLogger implements Loggable, Citable {

    public static final String TOTAL_COUNT_NAME = "totalChangeCount";
    public static final String COMPLETE_HISTORY_NAME = "completeHistory";

    public CompleteHistoryLogger(MarkovJumpsBeagleTreeLikelihood treeLikelihood) {
        this.tree = treeLikelihood.getTreeModel();
        this.patternCount = treeLikelihood.getPatternCount();

        treeTraitHistory = new TreeTrait[patternCount];
        for (int site = 0; site < patternCount; ++site) {
            String traitName = (patternCount == 0) ? MarkovJumpsBeagleTreeLikelihood.HISTORY : MarkovJumpsBeagleTreeLikelihood.HISTORY + "_" + (site + 1);
            treeTraitHistory[site] = treeLikelihood.getTreeTrait(traitName);
            if (treeTraitHistory[site] == null) {
                throw new RuntimeException("Tree '" + treeLikelihood.getId() + "' does not have a complete history trait at site " + (site + 1));
            }
        }

        treeTraitCount = treeLikelihood.getTreeTrait(MarkovJumpsBeagleTreeLikelihood.TOTAL_COUNTS);

        if (treeTraitCount == null) {
            throw new RuntimeException("No sum");
        }

        Logger.getLogger("dr.app.beagle").info("\tConstructing a complete history logger;  please cite:\n"
                + Citable.Utils.getCitationString(this));
    }

    public LogColumn[] getColumns() {

        LogColumn[] columns = new LogColumn[1 + patternCount];
        columns[0] = new LogColumn.Abstract(TOTAL_COUNT_NAME) {

            @Override
            protected String getFormattedValue() {
                return treeTraitCount.getTraitString(tree, null);
            }
        };

        for (int site = 0; site < patternCount; ++site) {

            String name = (patternCount == 0) ? COMPLETE_HISTORY_NAME : COMPLETE_HISTORY_NAME + "_" + (site + 1);
            final int anonSite = site;
            columns[1 + site] = new LogColumn.Abstract(name) {

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
                            String trait = treeTraitHistory[anonSite].getTraitString(tree, node);
                            if (trait.compareTo("{}") != 0) {
                                trait = trait.substring(1, trait.length() - 1);
                                StringTokenizer st = new StringTokenizer(trait, ",");
                                while (st.hasMoreTokens()) {
                                    if (!empty) {
                                        bf.append(",");
                                    }
                                    String event = st.nextToken();
                                    event = event.substring(1, event.length() - 1);
                                    StringTokenizer value = new StringTokenizer(event, ":");
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
        }
        return columns;
    }

    /**
     * @return a list of citations associated with this object
     */
    public List<Citation> getCitations() {
        List<Citation> citations = new ArrayList<Citation>();
        citations.add(CommonCitations.LEMEY_2012);
        citations.add(CommonCitations.SHAPIRO_2012);
        citations.add(CommonCitations.BLOOM_2012);
        return citations;
    }

    final private Tree tree;
    final private TreeTrait[] treeTraitHistory;
    final private TreeTrait treeTraitCount;
    final private int patternCount;
}
