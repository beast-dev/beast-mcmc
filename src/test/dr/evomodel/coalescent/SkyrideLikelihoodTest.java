package test.dr.evomodel.coalescent;

import dr.evolution.coalescent.TreeIntervalList;
import dr.evolution.io.Importer;
import dr.evolution.io.NewickImporter;
import dr.evolution.tree.Tree;
import dr.evomodel.coalescent.*;
import dr.evomodel.tree.DefaultTreeModel;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.Likelihood;
import dr.inference.model.Parameter;
import dr.math.MathUtils;
import junit.framework.TestCase;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * This is a very simplistic test case that ensues the skyride and skygrid likelihoods are unchanged by refactoring.
 * The values come from the master branch at commit id cf3d7370ca1a5b697f0f49be49765dcd6ad06dfb with class OldGRMRSkyride, GMRFMultilocusSkyride
 * which are assumed to be correct.
 */
public class SkyrideLikelihoodTest extends TestCase {
    public void setUp() throws Importer.ImportException, IOException {
        NewickImporter importer = new NewickImporter("(((0:0.5,(1:1.0,2:1.0)n6:1.0)n7:1.0,3:1.5)n8:1.0,(4:2.0,5:1.51)n9:1.5)n10;");
        MathUtils.setSeed(7);
        TreeModel tree = new DefaultTreeModel(importer.importTree(null));
        TreeIntervalList treeIntervals = new TreeIntervals(tree);
        List<TreeIntervalList> treeIntervalsList = new ArrayList<>();
        treeIntervalsList.add(treeIntervals);

//        List<Tree> treeList = new ArrayList<>();
//        treeList.add(tree);

        Parameter populationSizes = new Parameter.Default(1.0);
        Parameter groupSizes = new Parameter.Default(1.0);
        Parameter precision = new Parameter.Default(1.0);
        Parameter lambda = new Parameter.Default(1.0);

        Parameter populationGrid = new Parameter.Default(6,1.0);
        Parameter ploidyFactor = new Parameter.Default(1.0);

//        skyride = new OldGMRFSkyrideLikelihood(tree, populationSizes, groupSizes, precision, lambda, null, null, true, true);
//        skygrid = new GMRFMultilocusSkyrideLikelihood(treeList,populationGrid,groupSizes,precision,lambda,null,null,true,1.5,5,null,ploidyFactor);
//        gmrfSkyrideGradientNodes = new GMRFSkyrideGradient((OldGMRFSkyrideLikelihood)skyride,GMRFSkyrideGradient.WrtParameter.NODE_HEIGHTS, tree,null);
//        gmrfSkyrideGradientIntervals = new GMRFSkyrideGradient((OldGMRFSkyrideLikelihood)skyride,GMRFSkyrideGradient.WrtParameter.COALESCENT_INTERVAL, tree,null);

        skyride = new GMRFSkyrideLikelihood(treeIntervals, populationSizes, groupSizes, precision, lambda, null, null, true, true);
        skygrid = new GMRFMultilocusSkyrideLikelihood(treeIntervalsList,populationGrid,groupSizes,precision,lambda,null,null,true,1.5,5,null,ploidyFactor);

        gmrfSkyrideGradientNodes = new GMRFSkyrideGradient((GMRFSkyrideLikelihood)skyride,GMRFSkyrideGradient.WrtParameter.NODE_HEIGHTS, tree,null);
        gmrfSkyrideGradientIntervals = new GMRFSkyrideGradient((GMRFSkyrideLikelihood)skyride,GMRFSkyrideGradient.WrtParameter.COALESCENT_INTERVAL, tree,null);

    }

    public void testSkyride() {
        assertEquals(-13.837102559635337, skyride.getLogLikelihood(), 1e-10);
    }

    public void testSkygrid(){
        assertEquals(-14.756041059635336,skygrid.getLogLikelihood(),1e-10);
    }

    public void testGradientWRTnodes(){
        double[] gradLogDensity = {-1.103638323514327, -1.4715177646857693, -0.7357588823428847, -1.103638323514327, -0.36787944117144233};

        assertTrue(testArray(gradLogDensity,gmrfSkyrideGradientNodes.getGradientLogDensity(),1e-9));
    }
    public void testGradientWRTIntervals(){
        double[] gradLogDensity = {-4.782432735228751, -3.6787944117144233, -2.207276647028654, -1.103638323514327, -0.36787944117144233};
        assertTrue(testArray(gradLogDensity,gmrfSkyrideGradientIntervals.getGradientLogDensity(),1e-9));
    }

    private boolean testArray(double[] expected,double[] observed, double epsilon){
        if(expected.length!=observed.length){
            System.out.println(Arrays.toString(expected));
            System.out.println(Arrays.toString(observed));
            return false;
        }

        for(int i =0; i<expected.length;i++){
            if(observed[i]< expected[i]-epsilon|| observed[i]>expected[i]+epsilon){
                System.out.println(Arrays.toString(expected));
                System.out.println(Arrays.toString(observed));
                return false;
            }
        }
        return true;

    }


    private Likelihood skyride;
    private Likelihood skygrid;
    private GMRFSkyrideGradient gmrfSkyrideGradientNodes;
    private GMRFSkyrideGradient gmrfSkyrideGradientIntervals;
}
