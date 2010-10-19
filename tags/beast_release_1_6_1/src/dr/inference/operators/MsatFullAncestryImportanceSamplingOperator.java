package dr.inference.operators;

import dr.evomodel.tree.MicrosatelliteSamplerTreeModel;
import dr.evomodel.tree.TreeModel;
import dr.evomodel.substmodel.MicrosatelliteModel;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.inference.model.Parameter;
import dr.evolution.tree.Tree;
import dr.evolution.tree.NodeRef;
import dr.math.MathUtils;

/**
 * @author Chieh-Hsi Wu
 *
 * Produce an importance sample of the ancestry given a msat pattern and a tree.
 */
public class MsatFullAncestryImportanceSamplingOperator extends SimpleMCMCOperator{

    public static final String MSAT_FULL_ANCESTRY_IMPORTANCE_SAMPLING_OPERATOR = "MsatFullAncestryImportanceSamplingOperator";
    private Parameter parameter;
    private MicrosatelliteSamplerTreeModel msatSamplerTreeModel;
    private MicrosatelliteModel msatModel;
    private BranchRateModel branchRateModel;


    public MsatFullAncestryImportanceSamplingOperator(
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

        //get postOrder
        int[] postOrder = new int[tree.getNodeCount()];
        Tree.Utils.postOrderTraversalList(tree,postOrder);

        int extNodeCount = tree.getExternalNodeCount();
        double logq=0.0;
        for(int i = 0; i < postOrder.length; i ++){

            //if it's an internal node
            if(postOrder[i] >= extNodeCount){

                //getLikelihoodGiven the children
                NodeRef node = tree.getNode(postOrder[i]);
                NodeRef lc = tree.getChild(node,0);
                NodeRef rc = tree.getChild(node,1);
                int lcState = msatSamplerTreeModel.getNodeValue(lc);
                int rcState = msatSamplerTreeModel.getNodeValue(rc);
                double branchLeftLength = tree.getBranchLength(lc)*branchRateModel.getBranchRate(tree,lc);
                double branchRightLength = tree.getBranchLength(rc)*branchRateModel.getBranchRate(tree,rc);
                double[] probLbranch = msatModel.getColTransitionProbabilities(branchLeftLength, lcState);
                double[] probRbranch = msatModel.getColTransitionProbabilities(branchRightLength, rcState);
                double[] lik = new double[msatModel.getDataType().getStateCount()];
                int currState = (int)parameter.getParameterValue(msatSamplerTreeModel.getParameterIndexFromNodeNumber(postOrder[i]));
                //if node = root node
                if(i == postOrder.length -1){
                    //likelihood of root state also depends on the stationary distribution
                    double[] statDist = msatModel.getStationaryDistribution();
                    for(int j = 0; j < lik.length; j++){
                        lik[j] = probLbranch[j]*probRbranch[j]*statDist[j];
                    }

                }else{

                    for(int j = 0; j < lik.length; j++){
                        lik[j] = probLbranch[j]*probRbranch[j];
                    }

                }

                int sampledState = MathUtils.randomChoicePDF(lik);
                logq = logq + Math.log(lik[currState]) - Math.log(lik[sampledState]);
                parameter.setParameterValue(msatSamplerTreeModel.getParameterIndexFromNodeNumber(postOrder[i]),sampledState);
            }
        }
        
        return logq;
    }

    public String getPerformanceSuggestion(){
        return "None";
    }
    public String getOperatorName(){
        return MSAT_FULL_ANCESTRY_IMPORTANCE_SAMPLING_OPERATOR;
    }
    public int getStepCount(){
        return 1;
    }
}
