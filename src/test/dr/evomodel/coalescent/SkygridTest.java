package test.dr.evomodel.coalescent;

import dr.inference.model.MatrixParameter;
import dr.inference.model.Parameter;

import java.util.ArrayList;
import java.util.List;

import dr.evolution.coalescent.Coalescent;
import dr.evolution.coalescent.IntervalList;
import dr.evolution.coalescent.TreeIntervalList;
import dr.evolution.io.NewickImporter;
import dr.evolution.tree.Tree;
import dr.evolution.util.Units.Type;
import dr.evomodel.coalescent.CoalescentLikelihood;
import dr.evomodel.coalescent.GMRFMultilocusSkyrideLikelihood;
import dr.evomodel.coalescent.GMRFSkygridLikelihood;
import dr.evomodel.coalescent.GMRFSkyrideLikelihood;
import dr.evomodel.coalescent.TreeIntervals;
import dr.evomodel.coalescent.demographicmodel.ConstantPopulationModel;
import dr.evomodel.tree.DefaultTreeModel;
import dr.math.MathUtils;
import junit.framework.TestCase;

public class SkygridTest extends TestCase {

    Tree tree;
    List<Tree> treeList;
    TreeIntervals intervals;
    List<TreeIntervalList> intervalList;
    ConstantPopulationModel trueModel;
    Parameter popSize;
    CoalescentLikelihood likelihood;
    Parameter groupParameter = null;
    Parameter lambda = new Parameter.Default("LAMBDA", 1.0);
    Parameter phi = null;
    Parameter ploidyFactors= new Parameter.Default("PLOIDY", 1);;
    Parameter betaParameter= null;
    Parameter popParameter= null;
    MatrixParameter dMatrix = null;
    boolean timeAwareSmoothing = GMRFSkyrideLikelihood.TIME_AWARE_IS_ON_BY_DEFAULT;
    
    
    
    protected void setUp() throws Exception {
        NewickImporter importer = new NewickImporter("(((0:0.5,(1:1.0,2:1.0)n6:1.0)n7:1.0,3:1.5)n8:1.0,(4:2.0,5:1.51)n9:1.5)n10;");
        MathUtils.setSeed(7);
        
        tree = new DefaultTreeModel(importer.importTree(null));
        treeList = new ArrayList<Tree>();
        treeList.add(tree);

        intervals = new TreeIntervals(tree);

        intervalList = new ArrayList<TreeIntervalList>();
        intervalList.add(intervals);


        popSize = new Parameter.Default(1);
        popSize.setValue(0, 10.0);
        trueModel = new ConstantPopulationModel(popSize, Type.YEARS); 

        likelihood = new CoalescentLikelihood(intervals, trueModel);
    }

    public void testOldSkygrid(){
        Parameter cutOff = new Parameter.Default(1);
        cutOff.setValue(0,4.5 );
        Parameter numGridPoints = new Parameter.Default(2);
        Parameter precParameter = new Parameter.Default(1);
        precParameter.setValue(0, 0.1);

        popParameter = new Parameter.Default(2);
        double logPop = Math.log(popSize.getParameterValue(0));
        popParameter.setValue(0, logPop);
        popParameter.setValue(1, logPop);

        GMRFMultilocusSkyrideLikelihood skygrid = new GMRFMultilocusSkyrideLikelihood(treeList, popParameter, groupParameter, precParameter,
        lambda, betaParameter, dMatrix, timeAwareSmoothing, cutOff.getParameterValue(0), (int) numGridPoints.getParameterValue(0), phi, ploidyFactors);
        skygrid.getLogLikelihood();
        double logLikelihood = skygrid.peakLogCoalescentLikelihood();
        double constantLL  = likelihood.getLogLikelihood();
        assertEquals(constantLL,logLikelihood,1e-10);
    }
    
    public void testNewSkygrid(){
        Parameter cutOff = new Parameter.Default(1);
        cutOff.setValue(0,4.5 );
        Parameter numGridPoints = new Parameter.Default(2);
        Parameter precParameter = new Parameter.Default(1);
        precParameter.setValue(0, 0.1);

        popParameter = new Parameter.Default(2);
        double logPop = Math.log(popSize.getParameterValue(0));
        popParameter.setValue(0, logPop);
        popParameter.setValue(1, logPop);

        GMRFSkygridLikelihood skygrid = new GMRFSkygridLikelihood(intervalList, popParameter, groupParameter, precParameter,
        lambda, betaParameter, dMatrix, timeAwareSmoothing, cutOff.getParameterValue(0), (int) numGridPoints.getParameterValue(0), phi, ploidyFactors);

        skygrid.getLogLikelihood();
        double logLikelihood = skygrid.peakLogCoalescentLikelihood();
        double constantLL  = likelihood.getLogLikelihood();
        assertEquals(constantLL,logLikelihood,1e-10);
    }


    public void testOldSkygridRecentCutoff(){
        Parameter cutOff = new Parameter.Default(1);
        cutOff.setValue(0,2.0 );
        Parameter numGridPoints = new Parameter.Default(2);
        Parameter precParameter = new Parameter.Default(1);
        precParameter.setValue(0, 0.1);

        popParameter = new Parameter.Default(2);
        double logPop = Math.log(popSize.getParameterValue(0));
        popParameter.setValue(0, logPop);
        popParameter.setValue(1, logPop);

        GMRFMultilocusSkyrideLikelihood skygrid = new GMRFMultilocusSkyrideLikelihood(treeList, popParameter, groupParameter, precParameter,
        lambda, betaParameter, dMatrix, timeAwareSmoothing, cutOff.getParameterValue(0), (int) numGridPoints.getParameterValue(0), phi, ploidyFactors);
        skygrid.getLogLikelihood();
        double logLikelihood = skygrid.peakLogCoalescentLikelihood();
        double constantLL  = likelihood.getLogLikelihood();
        assertEquals(constantLL,logLikelihood,1e-10);
    }

    public void testNewSkygridRecentCutoff(){
        Parameter cutOff = new Parameter.Default(1);
        cutOff.setValue(0,2.0 );
        Parameter numGridPoints = new Parameter.Default(2);
        Parameter precParameter = new Parameter.Default(1);
        precParameter.setValue(0, 0.1);

        popParameter = new Parameter.Default(2);
        double logPop = Math.log(popSize.getParameterValue(0));
        popParameter.setValue(0, logPop);
        popParameter.setValue(1, logPop);

        GMRFSkygridLikelihood skygrid = new GMRFSkygridLikelihood(intervalList, popParameter, groupParameter, precParameter,
        lambda, betaParameter, dMatrix, timeAwareSmoothing, cutOff.getParameterValue(0), (int) numGridPoints.getParameterValue(0), phi, ploidyFactors);
        skygrid.getLogLikelihood();
        double logLikelihood = skygrid.peakLogCoalescentLikelihood();
        double constantLL  = likelihood.getLogLikelihood();
        assertEquals(constantLL,logLikelihood,1e-10);
    }



    public void testOldSkygridNoGridPoints(){
        Parameter cutOff = new Parameter.Default(1);
        cutOff.setValue(0,100.0 );
        Parameter numGridPoints = new Parameter.Default(2);
        Parameter precParameter = new Parameter.Default(1);
        precParameter.setValue(0, 0.1);

        popParameter = new Parameter.Default(2);
        double logPop = Math.log(popSize.getParameterValue(0));
        popParameter.setValue(0, logPop);
        popParameter.setValue(1, logPop);

        GMRFMultilocusSkyrideLikelihood skygrid = new GMRFMultilocusSkyrideLikelihood(treeList, popParameter, groupParameter, precParameter,
        lambda, betaParameter, dMatrix, timeAwareSmoothing, cutOff.getParameterValue(0), (int) numGridPoints.getParameterValue(0), phi, ploidyFactors);
        skygrid.getLogLikelihood();
        double logLikelihood = skygrid.peakLogCoalescentLikelihood();
        double constantLL  = likelihood.getLogLikelihood();
        assertEquals(constantLL,logLikelihood,1e-10);
    }



    public void testNewSkygridNoGridPoints(){
        Parameter cutOff = new Parameter.Default(1);
        cutOff.setValue(0,100.0 );
        Parameter numGridPoints = new Parameter.Default(2);
        Parameter precParameter = new Parameter.Default(1);
        precParameter.setValue(0, 0.1);

        popParameter = new Parameter.Default(2);
        double logPop = Math.log(popSize.getParameterValue(0));
        popParameter.setValue(0, logPop);
        popParameter.setValue(1, logPop);

        GMRFSkygridLikelihood skygrid = new GMRFSkygridLikelihood(intervalList, popParameter, groupParameter, precParameter,
        lambda, betaParameter, dMatrix, timeAwareSmoothing, cutOff.getParameterValue(0), (int) numGridPoints.getParameterValue(0), phi, ploidyFactors);
        skygrid.getLogLikelihood();
        double logLikelihood = skygrid.peakLogCoalescentLikelihood();
        double constantLL  = likelihood.getLogLikelihood();
        assertEquals(constantLL,logLikelihood,1e-10);
    }





    public void testOldSkygridMultifurcation(){
        Parameter cutOff = new Parameter.Default(1);
        cutOff.setValue(0,4.5 );
        Parameter numGridPoints = new Parameter.Default(2);
        Parameter precParameter = new Parameter.Default(1);
        precParameter.setValue(0, 0.1);

        popParameter = new Parameter.Default(2);
        double logPop = Math.log(popSize.getParameterValue(0));
        popParameter.setValue(0, logPop);
        popParameter.setValue(1, logPop);

        GMRFMultilocusSkyrideLikelihood skygrid = new GMRFMultilocusSkyrideLikelihood(treeList, popParameter, groupParameter, precParameter,
        lambda, betaParameter, dMatrix, timeAwareSmoothing, cutOff.getParameterValue(0), (int) numGridPoints.getParameterValue(0), phi, ploidyFactors);
        
        IntervalList intervals = skygrid.getTreeIntervals(0);

        ((dr.evolution.coalescent.TreeIntervals) intervals).setMultifurcationLimit(0.0);

        skygrid.getLogLikelihood();
        double logLikelihood = skygrid.peakLogCoalescentLikelihood();
        double constantLL  = likelihood.getLogLikelihood();

        assertEquals(constantLL,logLikelihood,1e-10);
        // Passing the test depends on multifurcation being -1.0
    }


}
