package dr.evomodel.speciation;

import dr.evolution.tree.NodeRef;
import dr.inference.distribution.ParametricDistributionModel;
import dr.inference.model.CompoundModel;
import dr.inference.model.Likelihood;
import dr.inference.model.Parameter;

/**
 * @author Joseph Heled
 *         Date: 8/12/2008
 */
public class SpeciesTreeBMPrior extends Likelihood.Abstract {

    private final SpeciesTreeModel sTree;
    //private final ParametricDistributionModel dist;
    private final ParametricDistributionModel tips;
    private final Parameter popSigma;
    private final Parameter stSigma;

    private static final double d1 = (1 - Math.exp(-1));
    private static final double f2 = -0.5 * Math.log(2*Math.PI);
    private final boolean logRoot;

//    public SpeciesTreeSimplePrior(SpeciesTreeModel sTree, ParametricDistributionModel dist, ParametricDistributionModel tipsPrior) {
//        super(new CompoundModel("STprior"));
//        this.sTree = sTree;
//        this.dist = dist;
//        this.tips = tipsPrior;
//
//        final CompoundModel cm = (CompoundModel)this.getModel();
//        cm.addModel(tipsPrior);
//        cm.addModel(dist);
//        cm.addModel(sTree);
//    }

    public SpeciesTreeBMPrior(SpeciesTreeModel sTree, Parameter popSigma, Parameter stSigma,
                              ParametricDistributionModel tipsPrior, boolean logRoot) {
        super(new CompoundModel("STBMprior"));

        this.sTree = sTree;
        this.popSigma = popSigma;
        this.stSigma = stSigma;
        this.tips = tipsPrior;
        this.logRoot = logRoot;

        final CompoundModel cm = (CompoundModel)this.getModel();

        cm.addModel(tipsPrior);
        cm.addModel(sTree);
    }

    protected double calculateLogLikelihood() {
        double logLike = 0;

        //if( true ) {
        final NodeRef root = sTree.getRoot();
        final double stRootHeight = sTree.getNodeHeight(root);

        final SpeciesTreeModel.RawPopulationHelper helper = sTree.getPopulationHelper();
        final double lim = stRootHeight ;//helper.geneTreesRootHeight() ;

        {
            final double sigma = stSigma.getParameterValue(0);
            final double s2 = 2 * sigma * sigma;
            int count = 0;
            double[] pops = new double[2];

            for(int nn = 0; nn < sTree.getNodeCount(); ++nn) {
                final NodeRef n = sTree.getNode(nn);

                if( sTree.isExternal(n) ) {
                    logLike += tips.logPdf(helper.tipPopulation(n));
                } else {
                    for(int nc = 0; nc < 2; ++nc) {
                        final NodeRef child = sTree.getChild(n, nc);
                        helper.getPopulations(n, nc, pops);
                        final double dt = sTree.getBranchLength(child) / lim;
                        final double pDiff = (pops[1] - pops[0]) / lim;
                        logLike -= pDiff * pDiff / (s2 * dt);
                        //need to adjest for dt!
                        logLike -= .5 * Math.log(dt);
                        count += 1;
                    }
                }
            }

            if( ! helper.perSpeciesPopulation() ) {
                helper.getRootPopulations(pops);
                final double dt = helper.geneTreesRootHeight()/lim - 1;
                //final double dt = (lim - stRootHeight) / lim;
                //  log(p1/lim) - log(p0/lim)  = log(p1/p0)
                final double pDiff = logRoot ? Math.log(pops[1]/pops[0]) : (pops[1] - pops[0])/lim;

                logLike -= pDiff * pDiff / (s2 * dt);
                //need to adjust for dt!
                logLike -= .5 * Math.log(dt);
                count += 1;
            }

            logLike += count * (f2 - Math.log(sigma));
        }

        if( helper.perSpeciesPopulation() ) {
            final double sigma = (popSigma != null ? popSigma : stSigma).getParameterValue(0);
            final double s2 = 2 * sigma * sigma;

            for(int ns = 0; ns < helper.nSpecies(); ++ns) {
                double[] times = helper.getTimes(ns);
                double[] pops = helper.getPops(ns);
                final double tMax = times[times.length-1];

                double ll = 0.0;
                double x = pops[0]/tMax;
                double t0 = 0.0;
                for(int k = 1; k < times.length; ++k) {
                    final double y = pops[k]/tMax;
                    final double dt = (times[k-1] - t0) / tMax;
                    ll -= (y - x)*(y-x) / (s2*dt);
                    // need to adjest for dt!
                    ll -= .5 * Math.log(dt);
                    x = y;
                    t0 = times[k-1];
                }
                ll += (times.length-1) * (f2 - Math.log(sigma));

                logLike += ll;
            }
        }
        //}
        return logLike;
    }

    protected boolean getLikelihoodKnown() {
		return false;
	}

}