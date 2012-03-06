package dr.evomodel.coalescent.operators;


import dr.evomodel.coalescent.GaussianProcessSkytrackLikelihood;
import dr.evomodel.coalescent.PopSizeStatistic;
import dr.evomodelxml.coalescent.operators.GaussianProcessSkytrackBlockUpdateOperatorParser;
import dr.inference.model.Parameter;
import dr.inference.operators.*;
import dr.math.MathUtils;
import no.uib.cipr.matrix.*;


import java.util.logging.Logger;

/* A Metropolis-Hastings operator to update the log population sizes and precision parameter jointly under a GP  prior
 *
 * @author Julia Palacios
 * @author Vladimir Minin
 * @author Marc Suchard
 * @version $Id: GMRFSkylineBlockUpdateOperator.java,v 1.5 2007/03/20 11:26:49 msuchard Exp $
 */

//TODO: create a new parameter "alpha" for the proposal on lambda
public class GaussianProcessSkytrackBlockUpdateOperator extends AbstractCoercableOperator {

    private double scaleFactor;
    private double lambdaScaleFactor;
//    private int fieldLength;

    private int maxIterations;
    private double stopValue;
    public static final double TWO_TIMES_PI =6.283185;
    private Parameter popSizeParameter;
    private Parameter precisionParameter;
    private Parameter lambdaParameter;
    private Parameter lambdaBoundParameter;
    private Parameter changePoints;
    private int [] GPcounts;
    private int [] GPtype;
    private double [] intervals;
    private double [] GPcoalfactor;
    private double [] coalfactor;
    private double alphaprior;
    private double betaprior;
    private SymmTridiagMatrix currentQ;
    private int numberPoints;  //Number of coalescent + latent points
    private double [] addPoints;
    private int [] add;
    private int fixedNumberPoints; //Number of coalescent or sampling points





    GaussianProcessSkytrackLikelihood GPvalue;
//    GMRFSkyrideLikelihood gmrfField;

    private double[] zeros;

    public GaussianProcessSkytrackBlockUpdateOperator(GaussianProcessSkytrackLikelihood GPLikelihood,
                                                      double weight, CoercionMode mode, double scaleFactor,
                                                      int maxIterations, double stopValue) {
        super(mode);
        GPvalue = GPLikelihood;     //before gmrfField
        popSizeParameter = GPLikelihood.getPopSizeParameter();
        changePoints=GPLikelihood.getChangePoints();
        GPcoalfactor=GPLikelihood.getGPcoalfactor();
        coalfactor=GPLikelihood.getcoalfactor();
        GPtype=GPLikelihood.getGPtype();
        GPcounts=GPLikelihood.getGPcounts();
        precisionParameter = GPLikelihood.getPrecisionParameter();
        lambdaParameter = GPLikelihood.getLambdaParameter();
        alphaprior=GPLikelihood.getAlphaParameter();
        betaprior=GPLikelihood.getBetaParameter();
        lambdaBoundParameter=GPLikelihood.getLambdaBoundParameter();
        currentQ=GPLikelihood.getWeightMatrix();
        numberPoints=popSizeParameter.getSize();
        fixedNumberPoints=GPcounts.length;
        int [] add = new int[fixedNumberPoints];
        double [] addPoints = new double[fixedNumberPoints];

        this.scaleFactor = scaleFactor;
        lambdaScaleFactor = 0.0;
//        fieldLength = popSizeParameter.getDimension();

        this.maxIterations = maxIterations;

        this.stopValue = stopValue;
        setWeight(weight);

        zeros = new double[numberPoints];

        System.err.println("At Operator start");


    }


    //change the 0.0005 to a parameter in BlockUpdate Parser
    private double getProposalLambda(double currentValue) {
        double proposal= MathUtils.uniform(currentValue-0.005,currentValue+0.005);
                //Symmetric proposal
        if (proposal<0){
         proposal=-proposal;
        }
        return proposal;
    }

    //change mixing values to parameters - see Eq. 7
    private double getPriorLambda(double lambdaMean, double epsilon, double lambdaValue){
        double res;
        if (lambdaValue < lambdaMean) {res=epsilon*(1/lambdaMean);}
        else {res=(1-epsilon)*(1/lambdaMean)*Math.exp(-(1/lambdaMean)*(lambdaValue-lambdaMean)); }
        return res;
    }


    private void  getNewUpperBound(double currentValue){
        double res;
        double newLambda, lambdaMean;
        newLambda = getProposalLambda(currentValue);
        lambdaMean = lambdaParameter.getParameterValue(0);
        double a=getPriorLambda(lambdaMean,0.01,newLambda)*(Math.pow(newLambda/currentValue,popSizeParameter.getSize()))*Math.exp((currentValue-newLambda)*GPvalue.getConstlik())/getPriorLambda(lambdaMean,0.01,currentValue);
        if (a>MathUtils.nextDouble()) {this.lambdaBoundParameter.setParameterValue(0,newLambda);
    }
    }

//    private int getNumberPoints(){
//        return numberPoints;
//    }

    private double getQuadraticForm(SymmTridiagMatrix currentQ, DenseVector x){
        DenseVector diagonal1 = new DenseVector(x.size());
        currentQ.mult(x, diagonal1);
        return x.dot(diagonal1);
    }
// Assumes the input vector x is ordered
    protected SymmTridiagMatrix getQmatrix(double precision, DenseVector x ) {
            SymmTridiagMatrix res;
            double trick=0.0;
            double[] offdiag = new double[x.size() - 1];
            double[] diag = new double[x.size()];



             for (int i = 0; i < x.size() - 1; i++) {
                    offdiag[i] = precision*(-1.0 / (x.get(i+1)-x.get(i)));
                 if (i< x.size()-2){
                    diag[i+1]= -offdiag[i]+precision*(1.0/(x.get(i+2)-x.get(i+1))+trick);
                 }
                }
//              Diffuse prior correction - intrinsic
             //Take care of the endpoints
            diag[0] = -offdiag[0]+precision*trick;

            diag[x.size() - 1] = -offdiag[x.size() - 2]+precision*(trick);
            res = new SymmTridiagMatrix(diag, offdiag);
        return res;
        }

    protected SymmTridiagMatrix getQmatrix(double precision, double[] x ) {
              SymmTridiagMatrix res;
              double trick=0.0;
              double[] offdiag = new double[x.length - 1];
              double[] diag = new double[x.length];



               for (int i = 0; i < x.length - 1; i++) {
                      offdiag[i] = precision*(-1.0 / (x[i+1]-x[i]));
                   if (i< x.length-2){
                      diag[i+1]= -offdiag[i]+precision*(1.0/(x[i+2]-x[i+1])+trick);
                   }
                  }
//              Diffuse prior correction - intrinsic
               //Take care of the endpoints
              diag[0] = -offdiag[0]+precision*trick;

              diag[x.length - 1] = -offdiag[x.length - 2]+precision*(trick);
              res = new SymmTridiagMatrix(diag, offdiag);
          return res;
          }






     protected QuadupleGP sortUpdate(double [] sortedData, double [] newData){
//         note that sortedData and newData are already ordered (minimum to maximum)
//         and last(sortedData) > last(newData)
          int newLength = sortedData.length + newData.length;
          double [] newList = new double [newLength];
          int [] newOrder = new int [newLength];

//          indexNew contains the index where the newData is stored (index ordered) in newList
//          indexOld contains the index where OldData is stored (index ordered) index newList
          int [] indexNew =new int[newData.length];
          int [] indexOld =new int[sortedData.length];

          int index2=sortedData.length;
          double pivot=newData[0];
          int k = 0;
          int k1 = 0;
          int k2=0;

          for (int j = 0; j < sortedData.length-1; j++){
             if (sortedData[j]<pivot) {
                 newOrder[k]=j;
                 newList[k]=sortedData[j];
                 indexOld[k2]=k;
                 k2+=1;
                 k+=1; }
             else {
                     if (index2<newLength){
                     newOrder[k]=index2;
                     index2+=1;
                     newList[k]=pivot;
                     indexNew[k1]=k;
                     pivot=newData[index2-sortedData.length];
                     k+=1;
                     k1+=1;}
                     else {
                         newOrder[k]=j;
                         newList[k]=sortedData[j];
                         indexOld[k2]=k;
                         k2+=1;
                         k+=1;
                     }

                 }
             }

      return new QuadupleGP(newList,newOrder,indexNew,indexOld);
     }


     protected Quaduple1GP sortUpdate(double [] sortedData, double newData){
//         note that sortedData and newData are already ordered (minimum to maximum)
//         and last(sortedData) > last(newData)
          int newLength = sortedData.length + 1;
          double [] newList = new double [newLength];
          int [] newOrder = new int [newLength];

//          indexNew contains the index where the newData is stored (index ordered) in newList
//          indexOld contains the index where OldData is stored (index ordered) index newList
          int  indexNew=0;
          int [] indexOld =new int[sortedData.length];

          int index2=sortedData.length;
          double pivot=newData;
          int k = 0;
          int k1 = 0;
          int k2=0;

          for (int j = 0; j < sortedData.length-1; j++){
             if (sortedData[j]<pivot) {
                 newOrder[k]=j;
                 newList[k]=sortedData[j];
                 indexOld[k2]=k;
                 k2+=1;
                 k+=1; }
             else {
                     if (index2<newLength){
                     newOrder[k]=index2;
                     index2+=1;
                     newList[k]=pivot;
                     indexNew=k;
                     pivot=newData;
                     k+=1;
                     k1+=1;}
                     else {
                         newOrder[k]=j;
                         newList[k]=sortedData[j];
                         indexOld[k2]=k;
                         k2+=1;
                         k+=1;
                     }

                 }
             }

      return new Quaduple1GP(newList,newOrder,indexNew,indexOld);
     }

    //Returns the index (position) of newData + Neighbors in the Ordered List
     protected int [] Neighbors(int [] indexNew, int numberTotalData){
      int [] Neighbors = new int[numberTotalData];
      int k=0;
      int latest=0;

         for (int j=0; j<indexNew.length;j++){

             if (indexNew[j]>latest) {
                 Neighbors[k]=indexNew[j]-1;
                 Neighbors[k+1]=indexNew[j];
                 Neighbors[k+2]=indexNew[j]+1;
                 k+=3;
             }
             if (indexNew[j]==latest & indexNew[j]>0){
                 Neighbors[k]=indexNew[j]+1;
                 k+=1;
             }
             if (indexNew[j]==0){Neighbors[k]=0;
                 Neighbors[k+1]=1;
                 k+=2;
             }
             latest=indexNew[j]+1;
          }
         int [] FinalNeighbors = new int[k+1];
         System.arraycopy(Neighbors,0,FinalNeighbors,0,k+1);
       return FinalNeighbors;
     }

      protected int [] Neighbors(int indexNew, int numberTotalData){
      int [] Neighbors = new int[3];
      int k=3;
             if (indexNew>0) {
                 Neighbors[0]=indexNew-1;
                 Neighbors[1]=indexNew;
                 Neighbors[2]=indexNew+1;

             }  else {
                 Neighbors[0]=0;
                 Neighbors[1]=1;
                 k=2;
             }


         int [] FinalNeighbors = new int[k];
         System.arraycopy(Neighbors,0,FinalNeighbors,0,k);
       return FinalNeighbors;
     }


//  Returns an array with the positions indicated in Index
    protected double [] SubsetData(double [] Data, int [] Index){
    double [] Res= new double [Index.length];
        for (int j=0;j<Index.length;j++){
            Res[j]=Data[Index[j]];
        }
        return Res;
  }

     protected double [] SubsetData(DenseVector Data, int [] Index){
     double [] Res= new double [Index.length];
        for (int j=0;j<Index.length;j++){
            Res[j]= Data.get(Index[j]);
        }
        return Res;
  }
    protected int [] SubsetData(int [] Data, int [] Index){
    int [] Res= new int [Index.length];
        for (int j=0;j<Index.length;j++){
            Res[j]=Data[Index[j]];
        }
        return Res;
  }

    protected PairIndex SubIndex (int [] Index, int numberOldData, int numberNewData){
        int [] newArray =new int[numberNewData];
        int [] oldArray =new int[numberOldData-numberNewData];
        int k=0;
        int k2=0;
        for (int j=0;j<Index.length;j++){
            if (Index[j]>=numberOldData){
                newArray[k]=Index[j];
                k+=1;
            } else {
                oldArray[k2]=Index[j];
                k2+=1;
            }
        }
     return new PairIndex(newArray,oldArray);
    }
//
//    protected QuadupleGP Neighbors(double [] sortedData, int [] orderOfData, int numberOldData){
////         sortedData is now the complete new data and orderOfData its order
////         numberOldData is the number of observations in the original data
////        Neighbors gives the subset of sortedData and orderOfData that contains newData and its neighbors
//        int newLength = sortedData.length;
//        int addLength=newLength-numberOldData;
//        int k=0;
//        int [] indicator = new int[newLength]; //indicator for neighbors
//        int sum=0;
//        int k1=0;
//        int k2=0;
//
//        if (orderOfData[0]>=numberOldData){indicator[0]=1; indicator[1]=0;}
//        for (int j =1; j<newLength-1;j++){
//
//          if (orderOfData[j]>=numberOldData){
//              indicator[j]=1;
//              indicator[j+1]=1;
//              indicator[j-1]=1;
//        }
//        }
//
//        for (int j=0;j<newLength;j++){
//         sum+=indicator[j];
//        }
//
//        double [] neighborData=new double[sum];
//        int [] neighborOrder=new int[sum];
//
//        int [] positionOld = new int[sum-addLength]; //indicator for OlData
//        int [] positionNew = new int[addLength]; //indicator for NewData
//
//
//        for (int j=0; j<newLength;j++){
//            if (indicator[j]==1){
//                neighborData[k]=sortedData[j];
//                neighborOrder[k]=orderOfData[j];
//                k+=1;
//                if (neighborOrder[k]>=numberOldData) {positionNew[k1]=k; k1+=1;} else {positionOld[k2]=k; k2+=1;}
//            }
//        }
//
//        return new QuadupleGP(neighborData,neighborOrder,positionNew, positionOld);
//    }


    public class PairGP{
        private double[] array1;
        private int[] array2;

        public PairGP (double [] array1, int [] array2){
            this.array1=array1;
            this.array2=array2;

        }
        public double[] getData() {return array1;}
        public int[] getOrder() {return array2;}
    }

      public class Pair1GP{
        private double array1;
        private int[] array2;

        public Pair1GP (double  array1, int [] array2){
            this.array1=array1;
            this.array2=array2;

        }
        public double getData() {return array1;}
        public int[] getOrder() {return array2;}
    }
       public class TripGP{
        private double [] array1;
        private int[] array2;
        private int[] array3;

        public TripGP (double []  array1, int [] array2, int[] array3){
            this.array1=array1;
            this.array2=array2;
            this.array3=array3;

        }
        public double [] getData() {return array1;}
        public int[] getOrder() {return array2;}
        public int[] getNewOrder() {return array3;}
    }

      public class Trip1GP{
        private double array1;
        private int[] array2;
        private int array3;

        public Trip1GP (double  array1, int [] array2, int array3){
            this.array1=array1;
            this.array2=array2;
            this.array3=array3;

        }
        public double getData() {return array1;}
        public int[] getOrder() {return array2;}
        public int getNewOrder() {return array3;}
    }

      public class PairIndex{
        private int[] array1;
        private int[] array2;

        public PairIndex (int [] array1, int [] array2){
            this.array1=array1;
            this.array2=array2;

        }
        public int[] getOrderNew() {return array1;}
        public int[] getOrderOld() {return array2;}
    }

    public class QuadupleGP{
        private double[] array1;
        private int[] array2;
        private int[] array3;
        private int[] array4;

        public QuadupleGP (double [] array1, int [] array2, int [] array3, int[] array4){
            this.array1=array1;
            this.array2=array2;
            this.array3=array3;
            this.array4=array4;

        }
        public double[] getData() {return array1;}
        public int[] getOrder() {return array2;}
        public int [] getPositionNew() {return array3;}
        public int [] getPositionOld() {return array4;}
    }

    public class Quaduple1GP{
        private double[] array1;
        private int[] array2;
        private int array3;
        private int[] array4;

        public Quaduple1GP (double [] array1, int [] array2, int array3, int[] array4){
            this.array1=array1;
            this.array2=array2;
            this.array3=array3;
            this.array4=array4;

        }
        public double[] getData() {return array1;}
        public int[] getOrder() {return array2;}
        public int  getPositionNew() {return array3;}
        public int [] getPositionOld() {return array4;}
    }
//    public class PairMatrices{
//        private double[] array1;
//        private int[] array2;
//
//        public PairMatrices (double [] array1, int [] array2){
//            this.array1=array1;
//            this.array2=array2;
//
//        }
//        public double[] getData() {return array1;}
//        public int[] getOrder() {return array2;}
//    }

    //Old function -delete

//    private double getNewLambda(double currentValue, double lambdaScale) {
//        double a = MathUtils.nextDouble() * lambdaScale - lambdaScale / 2;
//        double b = currentValue + a;
//        if (b > 1)
//            b = 2 - b;
//        if (b < 0)
//            b = -b;
//
//        return b;
//    }


     private double getNewPrecision(double currentValue, double quadraticTerm) {
//        double alphapost, betapost, alphaprior, betaprior;
        double alphaPost=  alphaprior+popSizeParameter.getSize()*0.5;
        double betaPost = betaprior+0.5*(1/currentValue)*quadraticTerm;
        return MathUtils.nextGamma(alphaPost,betaPost);
    }

//    private double getPrecision(double currentValue, double alpha, double beta) {
//            double length = scaleFactor - 1 / scaleFactor;
//            double returnValue;
//
//
//            if (scaleFactor == 1)
//                return currentValue;
//            if (MathUtils.nextDouble() < length / (length + 2 * Math.log(scaleFactor))) {
//                returnValue = (1 / scaleFactor + length * MathUtils.nextDouble()) * currentValue;
//            } else {
//                returnValue = Math.pow(scaleFactor, 2.0 * MathUtils.nextDouble() - 1) * currentValue;
//            }
//
//            return returnValue;
//        }


    //Old function -delete
//    private double getNewPrecision(double currentValue, double scaleFactor) {
//        double length = scaleFactor - 1 / scaleFactor;
//        double returnValue;
//
//
//        if (scaleFactor == 1)
//            return currentValue;
//        if (MathUtils.nextDouble() < length / (length + 2 * Math.log(scaleFactor))) {
//            returnValue = (1 / scaleFactor + length * MathUtils.nextDouble()) * currentValue;
//        } else {
//            returnValue = Math.pow(scaleFactor, 2.0 * MathUtils.nextDouble() - 1) * currentValue;
//        }
//
//        return returnValue;
//    }

    public DenseVector getMultiNormalMean(DenseVector CanonVector, UpperTriangBandMatrix CholeskyUpper) {

          DenseVector tempValue = new DenseVector(zeros);
          DenseVector Mean = new DenseVector(zeros);

//          UpperTriangBandMatrix CholeskyUpper = Cholesky.getU();

          // Assume Cholesky factorization of the precision matrix Q = LL^T

          // 1. Solve L\omega = b

          CholeskyUpper.transSolve(CanonVector, tempValue);

          // 2. Solve L^T \mu = \omega

          CholeskyUpper.solve(tempValue, Mean);

          return Mean;
      }


//    public DenseVector getMultiNormalMean(DenseVector CanonVector, BandCholesky Cholesky) {
//
//        DenseVector tempValue = new DenseVector(zeros);
//        DenseVector Mean = new DenseVector(zeros);
//
//        UpperTriangBandMatrix CholeskyUpper = Cholesky.getU();
//
//        // Assume Cholesky factorization of the precision matrix Q = LL^T
//
//        // 1. Solve L\omega = b
//
//        CholeskyUpper.transSolve(CanonVector, tempValue);
//
//        // 2. Solve L^T \mu = \omega
//
//        CholeskyUpper.solve(tempValue, Mean);
//
//        return Mean;
//    }


    public DenseVector getMultiNormal(DenseVector Mean, UpperTriangBandMatrix CholeskyUpper) {
        int length = Mean.size();
        DenseVector tempValue = new DenseVector(length);

        for (int i = 0; i < length; i++)
           tempValue.set(i, MathUtils.nextGaussian());

        DenseVector returnValue = new DenseVector(zeros);

//        UpperTriangBandMatrix CholeskyUpper = Cholesky.getU();
       // 3. Solve L^T v = z
        CholeskyUpper.solve(tempValue, returnValue);
        // 4. Return x = \mu + v
        returnValue.add(Mean);
        return returnValue;
    }

    public TripGP getGPvalues(double [] currentChangePoints, DenseVector currentGPvalues, double[] newChangePoints, double precision){
        int currentLength=currentChangePoints.length;
        int addLength=newChangePoints.length;
        int newLength = currentLength+addLength;
//       Takes the currentChangePoints and the newChangePoints into getData() but ordered and the order in getOrder()
//        assumes that order numbers 0..currentLength are the positions of currentChangePoints
//        and currentLength..newLength  - currentLentgth are the positions of newChangePoints
        QuadupleGP tempQuad= sortUpdate(currentChangePoints, newChangePoints);
//        index=tempPair.getOrder();
//        totalChangePoints=tempPair.getData();
//      This is a silly thing to only compute the Q.matrix for the neighbors of the newChangePoints
//      Takes the positions of the newData in the new complete sorted data and adds the positions of neighbors



//      Takes the positions of the newData in the new complete sorted data and adds the positions of neighbors
            int [] NeighborsIndex =Neighbors(tempQuad.getPositionNew(), newLength);
//      Retrieves the positions indicated in NeigborsIndex from the complete sorted data
            DenseVector tempData = new DenseVector(SubsetData(tempQuad.getData(),NeighborsIndex));
            SymmTridiagMatrix Q = getQmatrix(precision, tempData);
//        Retrieves the positions indicated in NeighborsIndex from the getOrder
            int [] NeighborsOriginal= SubsetData(tempQuad.getOrder(),NeighborsIndex);
//      Generates two arrays: one with the positions of newData and other with the positions of OldData in the Neighbors data TempData
            PairIndex Indicators = SubIndex(NeighborsOriginal,currentLength,addLength);



//        Matrix Qnew = Matrices.getSubMatrix(Q,tempNeighbor.getPositionNew(),tempNeighbor.getPositionNew());
//        Matrix Qother=Matrices.getSubMatrix(Q,tempNeighbor.getPositionNew(),tempNeighbor.getPositionOld());
            UpperSPDBandMatrix varf = new UpperSPDBandMatrix(Matrices.getSubMatrix(Q,Indicators.getOrderNew(),Indicators.getOrderNew()),1);
            BandCholesky U1 = new BandCholesky(addLength,1,true);
            U1.factor(varf);
            DenseVector part = new DenseVector(NeighborsOriginal.length-addLength);
            int [] GPpositions =SubsetData(NeighborsOriginal,Indicators.getOrderOld());
            DenseVector currentGPneighbors = new DenseVector(SubsetData(currentGPvalues,GPpositions));
            Matrices.getSubMatrix(Q,Indicators.getOrderNew(),Indicators.getOrderOld()).mult(currentGPneighbors,part);
            DenseVector mean = new DenseVector(getMultiNormalMean(part, U1.getU()));
            double [] addGPvalues = getMultiNormal(mean,U1.getU()).getData();
            return new TripGP(addGPvalues,tempQuad.getOrder(),tempQuad.getPositionNew());

    }


    public Trip1GP getGPvalues(double [] currentChangePoints, DenseVector currentGPvalues, double newChangePoints, double precision){
            int currentLength=currentChangePoints.length;
            int addLength=1;
            int newLength = currentLength+addLength;
//       Takes the currentChangePoints and the newChangePoints into getData() but ordered and the order in getOrder()
//        assumes that order numbers 0..currentLength are the positions of currentChangePoints
//        and currentLength..newLength  - currentLentgth are the positions of newChangePoints
            Quaduple1GP tempQuad= sortUpdate(currentChangePoints, newChangePoints);
//        index=tempPair.getOrder();
//        totalChangePoints=tempPair.getData();
//      This is a silly thing to only compute the Q.matrix for the neighbors of the newChangePoints
//      Takes the positions of the newData in the new complete sorted data and adds the positions of neighbors

        int [] NeighborsIndex =Neighbors(tempQuad.getPositionNew(), newLength);
//      Retrieves the positions indicated in NeigborsIndex from the complete sorted data
        DenseVector tempData = new DenseVector(SubsetData(tempQuad.getData(),NeighborsIndex));
        SymmTridiagMatrix Q = getQmatrix(precision, tempData);
//        Retrieves the positions indicated in NeighborsIndex from the getOrder
        int [] NeighborsOriginal= SubsetData(tempQuad.getOrder(),NeighborsIndex);
        double part=0.0;
        double varf=1.0;
        if (NeighborsOriginal.length==3){
           varf = Q.get(1,1);
           part = -currentGPvalues.get(NeighborsOriginal[0])*Q.get(1,0)-currentGPvalues.get(NeighborsOriginal[1])*Q.get(1,2);
        }
        if (NeighborsOriginal.length==2){
           varf = Q.get(0,0);
           part = -currentGPvalues.get(NeighborsOriginal[1])*Q.get(0,1);
        }
        double res =part/varf+MathUtils.nextGaussian()/Math.sqrt(varf);

        return new Trip1GP(res,tempQuad.getOrder(),tempQuad.getPositionNew());
        }

    public double sigmoidal (double x){
       return 1/(1+Math.exp(-x));
   }

//
    public void addOnePoint(double gpval, double location, int positions, int where, int [] index){
       GPcounts[where]+=1;
       popSizeParameter.addDimension(positions,gpval);
       changePoints.addDimension(positions,location);
       int [] GPtype2= new int[popSizeParameter.getSize()];
       double [] coalfactor2 = new double[popSizeParameter.getSize()];

        System.arraycopy(GPtype,0,GPtype2,0,positions);
        GPtype2[positions]=-1;
        System.arraycopy(GPtype,positions,GPtype2,positions+1,GPtype.length-positions);
        GPtype = (int []) GPtype2;


        System.arraycopy(coalfactor,0,coalfactor2,0,positions);
        coalfactor2[positions]=GPcoalfactor[where];
        System.arraycopy(coalfactor,positions,coalfactor2,positions+1,coalfactor.length-positions);
        coalfactor= (double []) coalfactor2;


    }

     public void delOnePoint(int positions, int where){
       GPcounts[where]-=1;
       popSizeParameter.removeDimension(positions);
       changePoints.removeDimension(positions);
       int [] GPtype2= new int[popSizeParameter.getSize()];
       double [] coalfactor2 = new double [popSizeParameter.getSize()];

       System.arraycopy(GPtype,0,GPtype2,0,positions);
       System.arraycopy(GPtype,positions+1,GPtype2,positions,GPtype2.length-positions);
       GPtype = (int []) GPtype2;

       System.arraycopy(coalfactor,0,coalfactor2,0,positions);
       System.arraycopy(coalfactor,positions+1,coalfactor2,positions,coalfactor2.length-positions);
       coalfactor= (double []) coalfactor2;
    }

    public void addMultiplePoints(double [] gpvals, double [] locations, int[] positions, int [] where, int [] index){
       int [] dummy = new int[where.length];
       double [] dummy2= new double[where.length];
       for (int j=0;j<where.length;j++){
        GPcounts[where[j]]+=1;
        popSizeParameter.addDimension(positions[j],gpvals[j]);
        changePoints.addDimension(positions[j], locations[j]);
        dummy[j]=-1;
        dummy2[j]=GPcoalfactor[where[j]];
       }
        int [] GPtype2 = new int [popSizeParameter.getSize()];
        double [] coalfactor2 = new double[popSizeParameter.getSize()];
        System.arraycopy(GPtype,0,GPtype2,0,GPtype.length);
        System.arraycopy(dummy,0,GPtype2,GPtype.length,GPtype2.length);
        System.arraycopy(coalfactor,0,coalfactor2,0,coalfactor.length);
        System.arraycopy(dummy2,0,coalfactor2,coalfactor.length,coalfactor2.length);
        GPtype=SubsetData(GPtype2, index);
        coalfactor=SubsetData(coalfactor2,index);

       }


    //Samples the number of thinned points and updates all the lists: GPType, changePoints, PopSize Parameter, coal.factor
    //and GPcounts. Proposes to add/delete a latent point in each intercoalescent interval
    public void numberThinned(double [] currentChangePoints, DenseVector currentPopSize, double currentPrecision) {
       double lowL = 0.0;
       double uppL =0.0;
       Trip1GP tempG;
       double accept=0.0;
       int where=0;
       int who=0;

//         Proposes to add/delete and proposes location uniformly in each intercoalescent interval
       for (int j=0; j<fixedNumberPoints; j++){

           if (MathUtils.nextDouble()<0.5) {
            uppL+=GPvalue.getInterval(j);
            addPoints[j]=MathUtils.uniform(lowL,uppL);
            lowL=uppL;
            tempG = getGPvalues(currentChangePoints, currentPopSize, addPoints[j],currentPrecision);
            accept=lambdaBoundParameter.getParameterValue(0)*GPvalue.getInterval(j)*sigmoidal(addPoints[j])*GPcoalfactor[j]/(GPcounts[j]+1);
               if (MathUtils.nextDouble()<accept) {
                   addOnePoint(tempG.getData(),addPoints[j],tempG.getNewOrder(),j,tempG.getOrder());
                   where+=1;
               }
           }
           else {
               if (GPcounts[j]>0) {
                who=where+MathUtils.nextInt(GPcounts[j]);
                accept=GPcounts[j]/(lambdaBoundParameter.getParameterValue(0)*GPvalue.getInterval(j)*sigmoidal(-popSizeParameter.getParameterValue(who))*coalfactor[who]);
                if (MathUtils.nextDouble()<accept){
                    delOnePoint(who,j);
                }
               }
           }

           where+=GPcounts[j];
       }
 }

    public void locationThinned(double [] currentChangePoints, DenseVector currentPopSize, double currentPrecision) {
//    Sample the interval
        double where;
        int who=0;
        double lowL = 0.0;
        double uppL =0.0;
        double cum=0.0;
        boolean got=false;
        double Tlen=0.0;
        for (int j=0; j<fixedNumberPoints; j++){
            if (GPcounts[j]>0) {
            Tlen+=GPvalue.getInterval(j);
            }
        }
//        Need to only consider the intevals with latent points
        int j=0;
        double addPts;
        int position = 0;
        double accept;
        where= MathUtils.uniform(0,Tlen);
        while (got==false) {
            position+=GPcounts[j];
            cum+=GPvalue.getInterval(j);
            if (GPcounts[j]>0){
                 uppL+=GPvalue.getInterval(j);
                 if (where>lowL & where<=uppL){
                     who=j;
                    got=true;
                    position-=GPcounts[j]-1;}
                 lowL=uppL;
            }
        j++;
        }
//Propose new points and get GPvalue
        Trip1GP tempG;
        for (int i=0; i<GPcounts[who]; i++){
         addPts=MathUtils.uniform(cum-GPvalue.getInterval(who),cum);
         tempG= getGPvalues(currentChangePoints, currentPopSize, addPts,currentPrecision);
         accept=sigmoidal(-tempG.getData())/sigmoidal(currentPopSize.get(position+i));
            if (MathUtils.nextDouble()<accept){
                addOnePoint(tempG.getData(),addPts,tempG.getNewOrder(),j,tempG.getOrder());
                delOnePoint(position+i,who);
            }
        }
    }

    public double forLikelihood(double [] Gvalues, int [] Gtype){
        double loglik=0.0;
        for (int j=0;j<Gvalues.length;j++){
          loglik-=Math.log(1+Math.exp(-Gtype[j]*Gvalues[j]));
        }
        return loglik;
    }

    public void sliceSampling(double [] currentChangePoints, DenseVector currentPopSize, double currentPrecision){
        double theta;
        double thetaPrime;
        int keep = 1;
        double [] zeross = new double[currentPopSize.size()];
        DenseVector v= new DenseVector(zeross);
        DenseVector v1 = new DenseVector(zeross);
        DenseVector v2 = new DenseVector(zeross);
        DenseVector proposal = new DenseVector(zeross);
        currentQ=getQmatrix(currentPrecision,currentChangePoints);

        UpperSPDBandMatrix U = new UpperSPDBandMatrix(currentQ,1);
        BandCholesky U1 = new BandCholesky(zeross.length,1,true);
        U1.factor(U);
        v= getMultiNormal(v1,U1.getU());
        theta=MathUtils.uniform(0,TWO_TIMES_PI);
        v1.add(Math.sin(theta),currentPopSize);
        v1.add(Math.cos(theta),v);
        v2.add(Math.cos(theta),currentPopSize);
        v2.add(-Math.sin(theta),v);

        double thetaMin=0.0;
        double thetaMax=TWO_TIMES_PI;
        double [] popSize = currentPopSize.getData();
        double loglik=Math.log(MathUtils.nextDouble())+forLikelihood(popSize,GPtype);
        double loglik2=0.0;
        while (keep==1){
            thetaPrime=MathUtils.uniform(thetaMin,thetaMax);
            proposal.add(Math.sin(thetaPrime),v1);
            proposal.add(Math.cos(thetaPrime),v2);
            double [] popSize2 = proposal.getData();
            loglik2-=forLikelihood(popSize2,GPtype);
            if (loglik2>loglik) {keep=2;}
            else {
                if (thetaPrime>theta) {thetaMin=thetaPrime;} else {thetaMax=thetaPrime;}
            }
        }
        for (int j=0; j<popSizeParameter.getSize();j++){
            popSizeParameter.setParameterValue(j,proposal.get(j));
        }

    }

//    public DenseVector getMultiNormal(DenseVector StandNorm, DenseVector Mean, BandCholesky Cholesky) {
//
//        DenseVector returnValue = new DenseVector(zeros);
//
//        UpperTriangBandMatrix CholeskyUpper = Cholesky.getU();
//
//        // 3. Solve L^T v = z
//
//        CholeskyUpper.solve(StandNorm, returnValue);
//
//        // 4. Return x = \mu + v
//
//        returnValue.add(Mean);
//
//        return returnValue;
//    }
//

//    public static DenseVector getMultiNormal(DenseVector Mean, UpperSPDDenseMatrix Variance) {
//        int length = Mean.size();
//        DenseVector tempValue = new DenseVector(length);
//        DenseVector returnValue = new DenseVector(length);
//        UpperSPDDenseMatrix ab = Variance.copy();
//
//        for (int i = 0; i < returnValue.size(); i++)
//            tempValue.set(i, MathUtils.nextGaussian());
//
//        DenseCholesky chol = new DenseCholesky(length, true);
//        chol.factor(ab);
//
//        UpperTriangDenseMatrix x = chol.getU();
//
//        x.transMult(tempValue, returnValue);
//        returnValue.add(Mean);
//        return returnValue;
//    }


    public static double logGeneralizedDeterminant(UpperTriangBandMatrix Matrix) {
        double returnValue = 0;

        for (int i = 0; i < Matrix.numColumns(); i++) {
            if (Matrix.get(i, i) > 0.0000001) {
                returnValue += Math.log(Matrix.get(i, i));
            }
        }

        return returnValue;
    }

//    public DenseVector newtonRaphson(double[] data, DenseVector currentGamma, SymmTridiagMatrix proposedQ) throws OperatorFailedException {
//        return newNewtonRaphson(data, currentGamma, proposedQ, maxIterations, stopValue);
//    }

//    public static DenseVector newNewtonRaphson(double[] data, DenseVector currentGamma, SymmTridiagMatrix proposedQ,
//                                               int maxIterations, double stopValue) throws OperatorFailedException {
//        DenseVector iterateGamma = currentGamma.copy();
//        DenseVector tempValue = currentGamma.copy();
//
//        int numberIterations = 0;
//
//
//        while (gradient(data, iterateGamma, proposedQ).norm(Vector.Norm.Two) > stopValue) {
//            try {
//                jacobian(data, iterateGamma, proposedQ).solve(gradient(data, iterateGamma, proposedQ), tempValue);
//            } catch (MatrixNotSPDException e) {
//                Logger.getLogger("dr.evomodel.coalescent.operators.GMRFSkyrideBlockUpdateOperator").fine("Newton-Raphson F");
//                throw new OperatorFailedException("");
//            } catch (MatrixSingularException e) {
//                Logger.getLogger("dr.evomodel.coalescent.operators.GMRFSkyrideBlockUpdateOperator").fine("Newton-Raphson F");
//
//                throw new OperatorFailedException("");
//            }
//            iterateGamma.add(tempValue);
//            numberIterations++;
//
//            if (numberIterations > maxIterations) {
//                Logger.getLogger("dr.evomodel.coalescent.operators.GMRFSkyrideBlockUpdateOperator").fine("Newton-Raphson F");
//                throw new OperatorFailedException("Newton Raphson algorithm did not converge within " + maxIterations + " step to a norm less than " + stopValue + "\n" +
//                        "Try starting BEAST with a more accurate initial tree.");
//            }
//        }
//
//        Logger.getLogger("dr.evomodel.coalescent.operators.GMRFSkyrideBlockUpdateOperator").fine("Newton-Raphson S");
//        return iterateGamma;
//
//    }

//    private static DenseVector gradient(double[] data, DenseVector value, SymmTridiagMatrix Q) {
//
//        DenseVector returnValue = new DenseVector(value.size());
//        Q.mult(value, returnValue);
//        for (int i = 0; i < value.size(); i++) {
//            returnValue.set(i, -returnValue.get(i) - 1 + data[i] * Math.exp(-value.get(i)));
//        }
//        return returnValue;
//    }
//
//
//    private static SPDTridiagMatrix jacobian(double[] data, DenseVector value, SymmTridiagMatrix Q) {
//        SPDTridiagMatrix jacobian = new SPDTridiagMatrix(Q, true);
//        for (int i = 0, n = value.size(); i < n; i++) {
//            jacobian.set(i, i, jacobian.get(i, i) + Math.exp(-value.get(i)) * data[i]);
//        }
//        return jacobian;
//    }


    //In GMRF there is only one MH block proposing precision + GMRF jointly. Here we have : numberThinned,
//    locationThinned, the GMRF (GP values), Gibbs sampling precision and the upper bound lambda
    public double doOperation() throws OperatorFailedException {

        System.err.println("Here doOperation starts");
//        System.exit(-1);

        double currentPrecision = this.precisionParameter.getParameterValue(0);
        DenseVector currentPopSize = new DenseVector(popSizeParameter.getParameterValues());
        double currentQuadratic = getQuadraticForm(currentQ, currentPopSize);
//       Gibbs sample new precision
        getNewPrecision(currentPrecision,currentQuadratic);
        currentPrecision=this.precisionParameter.getParameterValue(0);

//        proposes and updates lambdaBoundParameter
        double currentLambda = this.lambdaBoundParameter.getParameterValue(0);
        getNewUpperBound(currentLambda);
        currentLambda=this.lambdaBoundParameter.getParameterValue(0);

        double [] currentChangePoints = this.changePoints.getParameterValues();

//        numberThinned();



//        ArrayList<ComparableDouble> times = new ArrayList<ComparableDouble>();
//                ArrayList<Integer> childs = new ArrayList<Integer>();
//                collectAllTimes(tree, root, exclude, times, childs);
//                int[] indices = new int[times.size()];



//        double currentPrecision = precisionParameter.getParameterValue(0);
        double proposedPrecision = this.getNewPrecision(currentPrecision, scaleFactor);

//        double currentLambda = this.lambdaParameter.getParameterValue(0);
        double proposedLambda = currentLambda;

        precisionParameter.setParameterValue(0, proposedPrecision);
        lambdaParameter.setParameterValue(0, proposedLambda);

        DenseVector currentGamma = new DenseVector(GPvalue.getPopSizeParameter().getParameterValues());
        DenseVector proposedGamma;

//        SymmTridiagMatrix currentQ = GPvalue.getStoredScaledWeightMatrix(currentPrecision, currentLambda);
//        SymmTridiagMatrix proposedQ = GPvalue.getScaledWeightMatrix(proposedPrecision, proposedLambda);




//        double[] wNative = gmrfField.getSufficientStatistics();

//        UpperSPDBandMatrix forwardQW = new UpperSPDBandMatrix(proposedQ, 1);
//        UpperSPDBandMatrix backwardQW = new UpperSPDBandMatrix(currentQ, 1);
//
//        BandCholesky forwardCholesky = new BandCholesky(wNative.length, 1, true);
//        BandCholesky backwardCholesky = new BandCholesky(wNative.length, 1, true);
//
//        DenseVector diagonal1 = new DenseVector(fieldLength);
//        DenseVector diagonal2 = new DenseVector(fieldLength);
//        DenseVector diagonal3 = new DenseVector(fieldLength);
//
//        DenseVector modeForward = newtonRaphson(wNative, currentGamma, proposedQ.copy());
//
//        for (int i = 0; i < fieldLength; i++) {
//            diagonal1.set(i, wNative[i] * Math.exp(-modeForward.get(i)));
//            diagonal2.set(i, modeForward.get(i) + 1);
//
//            forwardQW.set(i, i, diagonal1.get(i) + forwardQW.get(i, i));
//            diagonal1.set(i, diagonal1.get(i) * diagonal2.get(i) - 1);
//        }
//
//        forwardCholesky.factor(forwardQW.copy());
//
//        DenseVector forwardMean = getMultiNormalMean(diagonal1, forwardCholesky);
//
//        DenseVector stand_norm = new DenseVector(zeros);
//
//        for (int i = 0; i < zeros.length; i++)
//            stand_norm.set(i, MathUtils.nextGaussian());
//
//        proposedGamma = getMultiNormal(stand_norm, forwardMean, forwardCholesky);
//
//
//        for (int i = 0; i < fieldLength; i++)
//            popSizeParameter.setParameterValueQuietly(i, proposedGamma.get(i));
//
//        ((Parameter.Abstract) popSizeParameter).fireParameterChangedEvent();
//

        double hRatio = 0;

//        diagonal1.zero();
//        diagonal2.zero();
//        diagonal3.zero();
//
//        DenseVector modeBackward = newtonRaphson(wNative, proposedGamma, currentQ.copy());
//
//        for (int i = 0; i < fieldLength; i++) {
//            diagonal1.set(i, wNative[i] * Math.exp(-modeBackward.get(i)));
//            diagonal2.set(i, modeBackward.get(i) + 1);
//
//            backwardQW.set(i, i, diagonal1.get(i) + backwardQW.get(i, i));
//            diagonal1.set(i, diagonal1.get(i) * diagonal2.get(i) - 1);
//        }
//
//        backwardCholesky.factor(backwardQW.copy());
//
//        DenseVector backwardMean = getMultiNormalMean(diagonal1, backwardCholesky);
//
//        for (int i = 0; i < fieldLength; i++) {
//            diagonal1.set(i, currentGamma.get(i) - backwardMean.get(i));
//        }
//
//        backwardQW.mult(diagonal1, diagonal3);
//
//        // Removed 0.5 * 2
//        hRatio += logGeneralizedDeterminant(backwardCholesky.getU()) - 0.5 * diagonal1.dot(diagonal3);
//        hRatio -= logGeneralizedDeterminant(forwardCholesky.getU() ) - 0.5 * stand_norm.dot(stand_norm);


//        return hRatio;
//        System.err.println("Prueba");
//        System.exit(-1);
        return 0;
    }

    //MCMCOperator INTERFACE
  // This is the only part where GPSBUOperateroParser is used
    public final String getOperatorName() {
        return GaussianProcessSkytrackBlockUpdateOperatorParser.BLOCK_UPDATE_OPERATOR;
    }

    public double getCoercableParameter() {
//        return Math.log(scaleFactor);
        return Math.sqrt(scaleFactor - 1);
    }

    public void setCoercableParameter(double value) {
//        scaleFactor = Math.exp(value);
        scaleFactor = 1 + value * value;
    }

    public double getRawParameter() {
        return scaleFactor;
    }

    public double getScaleFactor() {
        return scaleFactor;
    }

    public double getTargetAcceptanceProbability() {
        return 0.234;
    }

    public double getMinimumAcceptanceLevel() {
        return 0.1;
    }

    public double getMaximumAcceptanceLevel() {
        return 0.4;
    }

    public double getMinimumGoodAcceptanceLevel() {
        return 0.20;
    }

    public double getMaximumGoodAcceptanceLevel() {
        return 0.30;
    }

    public final String getPerformanceSuggestion() {

        double prob = MCMCOperator.Utils.getAcceptanceProbability(this);

//        double prob = Utils.getAcceptanceProbability(this);
        double targetProb = getTargetAcceptanceProbability();
        dr.util.NumberFormatter formatter = new dr.util.NumberFormatter(5);

        double sf = OperatorUtils.optimizeWindowSize(scaleFactor, prob, targetProb);

        if (prob < getMinimumGoodAcceptanceLevel()) {
            return "Try setting scaleFactor to about " + formatter.format(sf);
        } else if (prob > getMaximumGoodAcceptanceLevel()) {
            return "Try setting scaleFactor to about " + formatter.format(sf);
        } else return "";
    }


//  public DenseVector oldNewtonRaphson(double[] data, DenseVector currentGamma, SymmTridiagMatrix proposedQ) throws OperatorFailedException{
//  return newNewtonRaphson(data, currentGamma, proposedQ, maxIterations, stopValue);
//
//}
//
//public static DenseVector newtonRaphson(double[] data, DenseVector currentGamma, SymmTridiagMatrix proposedQ,
//int maxIterations, double stopValue) {
//
//DenseVector iterateGamma = currentGamma.copy();
//int numberIterations = 0;
//while (gradient(data, iterateGamma, proposedQ).norm(Vector.Norm.Two) > stopValue) {
//inverseJacobian(data, iterateGamma, proposedQ).multAdd(gradient(data, iterateGamma, proposedQ), iterateGamma);
//numberIterations++;
//}
//
//if (numberIterations > maxIterations)
//throw new RuntimeException("Newton Raphson algorithm did not converge within " + maxIterations + " step to a norm less than " + stopValue);
//
//return iterateGamma;
//}
//
//private static DenseMatrix inverseJacobian(double[] data, DenseVector value, SymmTridiagMatrix Q) {
//
//      SPDTridiagMatrix jacobian = new SPDTridiagMatrix(Q, true);
//      for (int i = 0; i < value.size(); i++) {
//          jacobian.set(i, i, jacobian.get(i, i) + Math.exp(-value.get(i)) * data[i]);
//      }
//
//      DenseMatrix inverseJacobian = Matrices.identity(jacobian.numRows());
//      jacobian.solve(Matrices.identity(value.size()), inverseJacobian);
//
//      return inverseJacobian;
//  }
}
