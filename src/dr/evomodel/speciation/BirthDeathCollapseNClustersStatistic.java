package dr.evomodel.speciation;

import dr.evolution.tree.NodeRef;
import dr.evomodelxml.speciation.BirthDeathCollapseNClustersStatisticParser;
import dr.inference.model.Statistic;

/**
 * @author Graham Jones
 *         Date: 01/10/2013
 */
public class BirthDeathCollapseNClustersStatistic extends Statistic.Abstract {
    private SpeciesTreeModel spptree;
    private BirthDeathCollapseModel bdcm;


    public BirthDeathCollapseNClustersStatistic(SpeciesTreeModel spptree, BirthDeathCollapseModel bdcm) {
        super(BirthDeathCollapseNClustersStatisticParser.BDC_NCLUSTERS_STATISTIC);
        this.spptree = spptree;
        this.bdcm = bdcm;
    }

    @Override
    public int getDimension() {
        return 1;
    }

    @Override
    public double getStatisticValue(int dim) {
        double eps = bdcm.getCollapseHeight();
        int ninodes = spptree.getInternalNodeCount();
        int n =  0;
        for (int i = 0; i < ninodes; i++) {
            double h = spptree.getNodeHeight(spptree.getInternalNode(i));
            if (h > eps) {
                n++;
            }
        }

        int nshort = 0;
        for (int i = 0; i < ninodes; i++) {
            NodeRef nr = spptree.getInternalNode(i);
            boolean collapse;
            if (spptree.isRoot(nr)) {
                collapse = (spptree.getNodeHeight(nr) < eps);
            } else {
                NodeRef anc  = spptree.getParent(nr);
                collapse = (spptree.getNodeHeight(nr) < eps  && spptree.getNodeHeight(anc) >= eps);
            }
            if (collapse) {
                nshort++;
            }
        }

        if (nshort  !=  n+1)
        {
            System.out.println("DEBUG in BirthDeathCollapseNClustersStatistic");
        }
        assert nshort  == n+1;

        return n+1;
    }
}
