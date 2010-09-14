package dr.inference.operators;

import dr.inference.model.Parameter;
import dr.evomodel.tree.MicrosatelliteSamplerTreeModel;
import dr.evomodel.tree.TreeModel;
import dr.evomodel.substmodel.MicrosatelliteModel;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evolution.tree.NodeRef;
import dr.math.MathUtils;

/**
 * @author Chieh-Hsi Wu
 *
 *
 */
public class MsatSingleAncestralStateGibbsOperator extends SimpleMCMCOperator implements GibbsOperator{
    public static final String MSAT_SINGLE_ANCESTAL_STATE_GIBBS_OPERATOR = "MsatSingleAncestralStateGibbsOperator";
    private Parameter parameter;
    private MicrosatelliteSamplerTreeModel msatSamplerTreeModel;
    private MicrosatelliteModel msatModel;
    private BranchRateModel branchRateModel;

    public MsatSingleAncestralStateGibbsOperator (
            Parameter parameter,
            MicrosatelliteSamplerTreeModel msatSamplerTreeModel,
            MicrosatelliteModel msatModel,
            BranchRateModel branchRateModel,
            double weight){

        super();
        this.parameter = parameter;
        this.msatSamplerTreeModel = msatSamplerTreeModel;
        this.msatModel = msatModel;
        this.branchRateModel = branchRateModel;
        setWeight(weight);
    }

    public double doOperation(){
        TreeModel tree = msatSamplerTreeModel.getTreeModel();

        
        int index = MathUtils.nextInt(parameter.getDimension());

        //double logq=0.0;


        //getLikelihoodGiven the children
        NodeRef node = tree.getNode(msatSamplerTreeModel.getNodeNumberFromParameterIndex(index));
        NodeRef lc = tree.getChild(node,0);
        NodeRef rc = tree.getChild(node,1);


        int lcState = msatSamplerTreeModel.getNodeValue(lc);
        int rcState = msatSamplerTreeModel.getNodeValue(rc);


        double branchLeftLength = tree.getBranchLength(lc)*branchRateModel.getBranchRate(tree,lc);
        double branchRightLength = tree.getBranchLength(rc)*branchRateModel.getBranchRate(tree,rc);


        double[] probLbranch = msatModel.getColTransitionProbabilities(branchLeftLength, lcState);
        double[] probRbranch = msatModel.getColTransitionProbabilities(branchRightLength, rcState);

        double[] lik = new double[msatModel.getDataType().getStateCount()];
        //int currState = (int)parameter.getParameterValue(index);
        //if node = root node

        if(tree.isRoot(node)){
            //likelihood of root state also depends on the stationary distribution
            double[] statDist = msatModel.getStationaryDistribution();
            for(int j = 0; j < lik.length; j++){
                lik[j] = probLbranch[j]*probRbranch[j]*statDist[j];
            }
        }else{
            NodeRef parent = tree.getParent(node);
            int pState = msatSamplerTreeModel.getNodeValue(parent);
            double branchParentLength = tree.getBranchLength(node)*branchRateModel.getBranchRate(tree,node);
            double[] probPbranch = msatModel.getRowTransitionProbabilities(branchParentLength,pState);
            for(int j = 0; j < lik.length; j++){
                lik[j] = probLbranch[j]*probRbranch[j]*probPbranch[j];
            }
        }

        int sampledState = MathUtils.randomChoicePDF(lik);
        //logq = logq + Math.log(lik[currState]) - Math.log(lik[sampledState]);
        parameter.setParameterValue(index,sampledState);
        return 0.0;
    }

    public String getPerformanceSuggestion(){
        return "None";
    }
    public String getOperatorName(){
        return MSAT_SINGLE_ANCESTAL_STATE_GIBBS_OPERATOR;
    }
    public int getStepCount(){
        return 1;
    }

}