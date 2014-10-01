package dr.app.beagle.evomodel.treelikelihood;
/**
 * ${CLASS_NAME}
 *
 * @author Andrew Rambaut
 * @version $Id$
 *
 * $HeadURL$
 *
 * $LastChangedBy$
 * $LastChangedDate$
 * $LastChangedRevision$
 */

import dr.evolution.alignment.PatternList;
import dr.evolution.datatype.DataType;
import dr.evomodel.tree.TreeModel;

public abstract class AbstractSinglePartitionTreeLikelihood extends AbstractTreeLikelihood {
    public AbstractSinglePartitionTreeLikelihood(String name, PatternList patternList, TreeModel treeModel) {
        super(name, treeModel);

        this.patternList = patternList;
        this.dataType = patternList.getDataType();
        patternCount = patternList.getPatternCount();
        stateCount = dataType.getStateCount();

        patternWeights = patternList.getPatternWeights();

    }

    /**
     * Set update flag for a pattern
     */
    protected void updatePattern(int i) {
        if (updatePattern != null) {
            updatePattern[i] = true;
        }
        likelihoodKnown = false;
    }

    /**
     * Set update flag for all patterns
     */
    protected void updateAllPatterns() {
        if (updatePattern != null) {
            for (int i = 0; i < patternCount; i++) {
                updatePattern[i] = true;
            }
        }
        likelihoodKnown = false;
    }

    public final double[] getPatternWeights() {
        return patternWeights;
    }

    /**
     * Forces a complete recalculation of the likelihood next time getLikelihood is called
     */
    public void makeDirty() {
        super.makeDirty();
        updateAllPatterns();
    }

    /**
     * the patternList
     */
    protected PatternList patternList = null;
    protected DataType dataType = null;

    /**
     * the pattern weights
     */
    protected double[] patternWeights;

    /**
     * the number of patterns
     */
    protected int patternCount;

    /**
     * the number of states in the data
     */
    protected int stateCount;

    /**
     * Flags to specify which patterns are to be updated
     */
    protected boolean[] updatePattern = null;


}
