package dr.evomodel.approxPopTree;

import dr.evolution.alignment.Patterns;
import dr.evolution.tree.NodeRef;
import dr.inference.model.AbstractModel;

import java.util.LinkedList;

/**
 * Package: AbstractPopulationMRCAModel
 * Description:
 * <p/>
 * <p/>
 * Created by
 *
 * @author Alexander V. Alekseyenko (alexander.alekseyenko@gmail.com)
 *         Date: Feb 1, 2010
 *         Time: 5:30:01 PM
 */
abstract public class AbstractPopulationMRCAModel extends AbstractModel {
    protected double time; // time of the population node above
    protected double tMRCA;

    public AbstractPopulationMRCAModel(String name, double populationTime) {
        super(name);
        time = populationTime;
    }

    public double getMRCATime(LinkedList<NodeRef> nodes) { // get time of actual MRCA for the sequence nodes below
        return tMRCA;
    }

    abstract public double drawMRCATime(LinkedList<NodeRef> nodes); // re-draw from the distribution of t(MRCA)

    abstract public double[][] getMRCAPartials(LinkedList<NodeRef> nodes, Patterns patterns); // get partials at Population node

    abstract public double[][] drawMRCAPartials(LinkedList<NodeRef> nodes, Patterns patterns); // re-draw partial at Population node

    protected double[][] computeProfilePartials(LinkedList<NodeRef> nodes, Patterns patterns) { // compute profile partials for the sequence nodes, useful for the above methods

        int patternCount = patterns.getPatternCount();
        int stateCount = patterns.getStateCount();
        double partials[][] = new double[patternCount][stateCount];
        // Do the profile computation here
        //
        return partials;
    }
}
