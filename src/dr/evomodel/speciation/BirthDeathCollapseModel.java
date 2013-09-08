package dr.evomodel.speciation;

/**
 * Created with IntelliJ IDEA.
 * User: Graham
 * Date: 01/09/13
 * Time: 14:27
 * To change this template use File | Settings | File Templates.
 */


import dr.evolution.util.Taxon;
import dr.evolution.util.Units;
import dr.evolution.tree.Tree;
import dr.inference.model.Parameter;

import java.util.Set;



public class BirthDeathCollapseModel extends SpeciationModel {
    private Parameter birthDiffRate; // lambda - mu
    private Parameter relativeDeathRate; // mu/lambda
    private Parameter originHeight;
    private final double collapseWeight;
    private final double collapseHeight;

    public BirthDeathCollapseModel(String modelName, Tree tree, Parameter birthDiffRate, Parameter relativeDeathRate,
                                   Parameter originHeight, double collW, double collH, Units.Type units) {
        super(modelName, units);
        this.collapseWeight = collW;
        this.collapseHeight = collH;

        this.birthDiffRate = birthDiffRate;
        addVariable(birthDiffRate);
        birthDiffRate.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, 1));

        this.relativeDeathRate = relativeDeathRate;
        addVariable(relativeDeathRate);
        relativeDeathRate.addBounds(new Parameter.DefaultBounds(1.0, 0.0, 1));

        this.originHeight = originHeight;
        originHeight.setParameterValue(0, 1.05 * tree.getNodeHeight(tree.getRoot()));
        addVariable(originHeight);
        originHeight.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, 1));
    }




    @Override
    public double calculateTreeLogLikelihood(Tree tree) {
        double logpt = 0.0;
        int ninodes = tree.getInternalNodeCount();
        int ntips = tree.getExternalNodeCount();
        double alpha = birthDiffRate.getParameterValue(0);
        double beta = relativeDeathRate.getParameterValue(0);
        double tor = originHeight.getParameterValue(0);
        double w = collapseWeight;

        double rooth = tree.getNodeHeight(tree.getRoot());
        if (rooth > tor) {
            return Double.NEGATIVE_INFINITY;
        }

        logpt += originHeightLogLikelihood(tor, alpha, beta, w, ntips);

        for (int n = 0; n < ninodes; n++) {
            final double height = tree.getNodeHeight(tree.getInternalNode(n));
            double usualpn = nodeHeightLikelihood(height, tor, alpha, beta);
            double collapsedpn = (height < collapseHeight) ? 1.0 / collapseHeight : 0.0;
            logpt += Math.log((1.0 - w) * usualpn + w * collapsedpn);
        }

        return logpt;
    }


    private double originHeightLogLikelihood(double t, double a, double b, double w, int n) {
        double E = Math.exp(-a * t);
        double B = (1 - E) / (1-b*E);
        double z = 0.0;
        z += Math.log(a);
        z += Math.log(1 - b);
        z -= a * t;
        z -= 2 * Math.log(1 - b * E);
        z += (n-2) * Math.log(w + (1 - w) * B);
        z +=  Math.log(w + n * (1 - w) * B);
        return z;
    }


    private double nodeHeightLikelihood(double s, double t, double a, double b) {
        double Es = Math.exp(-a * s);
        double Et = Math.exp(-a * t);
        double z = 0.0;
        if (s < t) {
            z = a;
            z *= (1 - b);
            z *= Es;
            z /= (1 - b * Es) * (1 - b * Es);
            z *= (1 - b * Et);
            z /= (1 - Et);
        }
        return z;
    }


    @Override
    public double calculateTreeLogLikelihood(Tree tree, Set<Taxon> exclude) {
        // not implemented.
        return Double.NEGATIVE_INFINITY;
    }
}
