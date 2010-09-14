package dr.evomodel.MSSD;

import dr.evolution.alignment.PatternList;
import dr.evolution.util.Taxon;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.sitemodel.SiteModel;
import dr.evomodel.tree.TreeModel;
import dr.evomodelxml.MSSD.SingleTipObservationProcessParser;
import dr.inference.model.Parameter;

/**
 * Package: SingleTipObservationProcess
 * Description:
 * <p/>
 * <p/>
 * Created by
 * Alexander V. Alekseyenko (alexander.alekseyenko@gmail.com)
 * Date: Feb 19, 2008
 * Time: 2:57:14 PM
 */
public class SingleTipObservationProcess extends AnyTipObservationProcess {
    protected Taxon sourceTaxon;

    public SingleTipObservationProcess(TreeModel treeModel, PatternList patterns, SiteModel siteModel,
                                       BranchRateModel branchRateModel, Parameter mu, Parameter lam, Taxon sourceTaxon) {
        super(SingleTipObservationProcessParser.MODEL_NAME, treeModel, patterns, siteModel, branchRateModel, mu, lam);
        this.sourceTaxon = sourceTaxon;
    }

    public double calculateLogTreeWeight() {
        return -lam.getParameterValue(0) / (getAverageRate() * mu.getParameterValue(0));
    }

}
