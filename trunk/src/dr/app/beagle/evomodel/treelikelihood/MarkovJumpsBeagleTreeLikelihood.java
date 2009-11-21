package dr.app.beagle.evomodel.treelikelihood;

import dr.app.beagle.evomodel.substmodel.SubstitutionModel;
import dr.app.beagle.evomodel.substmodel.MarkovJumpsSubstitutionModel;
import dr.app.beagle.evomodel.sitemodel.SiteRateModel;
import dr.app.beagle.evomodel.sitemodel.BranchSiteModel;
import dr.evolution.datatype.DataType;
import dr.evolution.alignment.PatternList;
import dr.inference.model.MatrixParameter;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.tree.TreeModel;

/**
 * @author Marc Suchard
 * @author Vladimir Minin
 *         <p/>
 *         A base class for implementing Markov chain-induced counting processes (markovjumps) in BEAST using BEAGLE
 *         This work is supported by NSF grant 0856099
 *         <p/>
 *         Minin VN and Suchard MA (2008) Counting labeled transitions in continous-time Markov models of evolution.
 *         Journal of Mathematical Biology, 56, 391-412.
 */
public class MarkovJumpsBeagleTreeLikelihood extends AncestralStateBeagleTreeLikelihood {


    public MarkovJumpsBeagleTreeLikelihood(PatternList patternList, TreeModel treeModel,
                                           BranchSiteModel branchSiteModel, SiteRateModel siteRateModel,
                                           BranchRateModel branchRateModel, boolean useAmbiguities,
                                           PartialsRescalingScheme scalingScheme, DataType dataType, String tag,
                                           SubstitutionModel substModel, MatrixParameter registerMatrixParameter) {
        super(patternList, treeModel, branchSiteModel, siteRateModel, branchRateModel, useAmbiguities,
                scalingScheme, dataType, tag, substModel);

        markovjumps = new MarkovJumpsSubstitutionModel(substModel);
    }

    private MarkovJumpsSubstitutionModel markovjumps;
}
