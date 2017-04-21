/*
 * GaussianProcessSkytrackBlockUpdateOperator.java
 *
 * Copyright (c) 2002-2015 Alexei Drummond, Andrew Rambaut and Marc Suchard
 *
 * This file is part of BEAST.
 * See the NOTICE file distributed with this work for additional
 * information regarding copyright ownership and licensing.
 *
 * BEAST is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 *  BEAST is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with BEAST; if not, write to the
 * Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
 * Boston, MA  02110-1301  USA
 */

package dr.evomodel.coalescent.operators;


//import com.sun.servicetag.SystemEnvironment;
import dr.evomodel.coalescent.GaussianProcessSkytrackLikelihood;

import dr.evomodelxml.coalescent.operators.GaussianProcessSkytrackBlockUpdateOperatorParser;
//import dr.inference.mcmc.MCMCOptions;
import dr.inference.model.Parameter;
import dr.inference.operators.*;
import dr.math.MathUtils;
//import dr.math.matrixAlgebra.SymmetricMatrix;
import no.uib.cipr.matrix.*;

import java.io.BufferedWriter;
//import java.io.FileNotFoundException;
import java.io.FileWriter;

/* A Metropolis-Hastings/Gibbs operator to update the log population sizes and precision parameter jointly under a GP  prior
 *
 * @author Julia Palacios
 * @author Vladimir Minin
 * @author Marc Suchard
 * @version $Id: GMRFSkylineBlockUpdateOperator.java,v 1.5 2007/03/20 11:26:49 msuchard Exp $
 */

//public class GaussianProcessSkytrackBlockUpdateOperator extends AbstractCoercableOperator{

public class GaussianProcessSkytrackBlockUpdateOperator extends SimpleMCMCOperator implements GibbsOperator{
    private double scaleFactor;
//    private double lambdaScaleFactor;
//    private int fieldLength;

//    private int maxIterations;
//    private double stopValue;
    public static final double TWO_TIMES_PI =6.283185;
    private Parameter popSizeParameter;
    private Parameter precisionParameter;
    private Parameter lambdaParameter;
    private Parameter lambdaBoundParameter;
    private Parameter changePoints;
    private Parameter GPcounts;
    private Parameter GPtype;
//    private Parameter popValue;   //This will actually have the Effective Pop. Size for the grid
    private Parameter CoalCounts;
    private Parameter numPoints;
//    private double [] intervals;
    private double [] GPcoalfactor;
    private Parameter coalfactor;
    private double alphaprior;
    private double betaprior;
    private SymmTridiagMatrix currentQ;
    private int numberPoints;  //Number of coalescent  points
//    private double [] addPoints;
//    private int [] add;
//    private int fixedNumberPoints; //Number of coalescent or sampling points
    private int [] CoalPosIndicator;
    private double [] CoalTime;



    GaussianProcessSkytrackLikelihood GPvalue;
//    GMRFSkyrideLikelihood gmrfField;

//    private double[] zeros;

    public GaussianProcessSkytrackBlockUpdateOperator(GaussianProcessSkytrackLikelihood GPLikelihood,
                                                      double weight, CoercionMode mode, double scaleFactor,
                                                      int maxIterations, double stopValue) {
//        super(mode);
        GPvalue = GPLikelihood;     //before gmrfField
        popSizeParameter = GPLikelihood.getPopSizeParameter();

//        popValue=GPLikelihood.getPopValue();
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

        CoalPosIndicator=GPLikelihood.getCoalPosIndicator();
        CoalTime=GPLikelihood.getCoalTime();
//        fixedNumberPoints=GPcounts.getSize();
        CoalCounts=GPLikelihood.getCoalCounts();
        numberPoints=CoalCounts.getSize();
        numPoints=GPLikelihood.getNumPoints();
//        int [] add = new int[fixedNumberPoints];
//        double [] addPoints = new double[fixedNumberPoints];

        this.scaleFactor = scaleFactor;
//        lambdaScaleFactor = 0.0;
//        fieldLength = popSizeParameter.getDimension();

//        this.maxIterations = maxIterations;

//        this.stopValue = stopValue;
        setWeight(weight);

//        zeros = new double[numberPoints];





    }


    public GaussianProcessSkytrackBlockUpdateOperator() {

    }
    //change the 0.0005 to a parameter in BlockUpdate Parser
    private double getProposalLambda(double currentValue) {
        double proposal= MathUtils.uniform(currentValue-0.00001,currentValue+0.00001);
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
            double trick=0;
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
              double trick=0.00000000001;
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
          double pivot2=newData[0];
          double pivot1=sortedData[0];
          int k1 = 0;

          for (int j = 0; j < newLength; j++){
             if (index2<newLength) {
                 if (pivot1<pivot2) {
                     newList[j]=pivot1;
                     newOrder[j]=k1;
                     indexOld[k1]=j;
                     k1++;
                     pivot1=sortedData[k1];
                 } else {
                     newList[j]=pivot2;

                     newOrder[j]=index2;

                     indexNew[index2-sortedData.length]=j;

                     index2++;
                     if (index2<newLength){
                         pivot2=newData[index2-sortedData.length];
                     }
                 }

             } else {
                 newList[j]=sortedData[k1];
                 newOrder[j]=k1;
                 indexOld[k1]=j;
                 k1++;
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
          double pivot2=newData;
          double pivot1=sortedData[0];
          int k1 = 0;

         for (int j = 0; j < newLength; j++){
             if (index2<newLength) {
                 if (pivot1<pivot2) {
                     newList[j]=pivot1;
                     newOrder[j]=k1;
                     indexOld[k1]=j;
                     k1++;
                     pivot1=sortedData[k1];
                 } else {
                     newList[j]=pivot2;

                     newOrder[j]=index2;

                     indexNew=j;

                     index2++;
                     if (index2<newLength){
                         pivot2=newData;
                     }
                 }

             } else {
                 newList[j]=sortedData[k1];
                 newOrder[j]=k1;
                 indexOld[k1]=j;
                 k1++;
             }
         }


      return new Quaduple1GP(newList,newOrder,indexNew,indexOld);
     }

    //Returns the index (position) of newData + Neighbors in the Ordered List
     protected int [] Neighbors(int [] indexNew, int numberTotalData){
      int [] Neighbors = new int[numberTotalData];
      int k=0;

         for (int j=0; j<indexNew.length-1;j++){

             if (indexNew[j]==0) {
                 if (indexNew[j+1]>2){
                     Neighbors[k]=0;
                     Neighbors[k+1]=1;
                     k+=2;
                 }
                 if (indexNew[j+1]==2){
                     Neighbors[k]=0;
                     k++;
                 }
             } else {
                 if ((indexNew[j+1]-indexNew[j])>2){
                     Neighbors[k]=indexNew[j]-1;
                     Neighbors[k+1]=indexNew[j];
                     Neighbors[k+2]=indexNew[j]+1;
                     k+=3;
                 }
                 if ((indexNew[j+1]-indexNew[j])==2){
                     Neighbors[k]=indexNew[j]-1;
                     Neighbors[k+1]=indexNew[j];
                     k+=2;
                 }
                 if ((indexNew[j+1]-indexNew[j])==1){
                     Neighbors[k]=indexNew[j]-1;
                     k++;
                 }

             }
         }
         Neighbors[k]=indexNew[indexNew.length-1]-1;
         Neighbors[k+1]=indexNew[indexNew.length-1];
         Neighbors[k+2]=indexNew[indexNew.length-1]+1;
         k+=3;

         int [] FinalNeighbors = new int[k];
         System.arraycopy(Neighbors,0,FinalNeighbors,0,k);
       return FinalNeighbors;
     }

      protected int [] Neighbors(int indexNew){
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

//    Returns the index in the Neighbors data that contain new data
    protected PairIndex SubIndex (int [] Index, int numberOldData, int numberNewData){
        int [] newArray =new int[numberNewData];
        int [] oldArray =new int[numberOldData];
        int k=0;
        int k2=0;
        for (int j=0;j<Index.length;j++){
            if (Index[j]>=numberOldData){
                newArray[k]=j;
                k+=1;
            } else {
                oldArray[k2]=j;
                k2+=1;
            }
        }
        int [] oldArray2 = new int[k2];
        System.arraycopy(oldArray,0,oldArray2,0,k2);
     return new PairIndex(newArray,oldArray2);
    }

    private class PairGP{
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

      private class Trip1GP{
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

     private double getNewPrecision(double currentValue, double quadraticTerm) {

         double alphaPost=  alphaprior+popSizeParameter.getSize()*0.5;

         double betaPost = betaprior+0.5*(1/currentValue)*quadraticTerm;

         return MathUtils.nextGamma(alphaPost,betaPost);
    }

    public DenseVector getMultiNormalMean(DenseVector CanonVector, UpperTriangBandMatrix CholeskyUpper) {

          DenseVector tempValue = new DenseVector(CanonVector.size());
          DenseVector Mean = new DenseVector(CanonVector.size());

//            UpperTriangBandMatrix CholeskyUpper = Cholesky.getU();

          // Assume Cholesky factorization of the precision matrix Q = LL^T

          // 1. Solve L\omega = b
//        System.err.println(CanonVector.size()+"and"+tempValue.size());


          CholeskyUpper.transSolve(CanonVector, tempValue);


//

        // 2. Solve L^T \mu = \omega

          CholeskyUpper.solve(tempValue, Mean);

          return Mean;
      }




    public DenseVector getMultiNormal(DenseVector Mean, UpperTriangBandMatrix CholeskyUpper) {
        int length = Mean.size();
        DenseVector tempValue = new DenseVector(length);

        for (int i = 0; i < length; i++)
           tempValue.set(i, MathUtils.nextGaussian());

        DenseVector returnValue = new DenseVector(Mean.size());

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
//          index=tempPair.getOrder();
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

            UpperSPDBandMatrix varf = new UpperSPDBandMatrix(Matrices.getSubMatrix(Q,Indicators.getOrderNew(),Indicators.getOrderNew()),1);
            BandCholesky U1 = new BandCholesky(addLength,1,true);

            U1.factor(varf);


            DenseVector part = new DenseVector(addLength);
            int [] GPpositions =SubsetData(NeighborsOriginal,Indicators.getOrderOld());
            DenseVector currentGPneighbors = new DenseVector(SubsetData(currentGPvalues,GPpositions));

            Matrix first = Matrices.getSubMatrix(Q,Indicators.getOrderNew(),Indicators.getOrderOld());

            first.mult(-1,currentGPneighbors,part);

            DenseVector mean = new DenseVector(getMultiNormalMean(part, U1.getU()));

        double [] addGPvalues = getMultiNormal(mean,U1.getU()).getData();


        return new TripGP(addGPvalues,tempQuad.getOrder(),tempQuad.getPositionNew());

    }


    public double[] getGPvaluesS(double[] currentChangePoints, double[] currentGPvalues1, double[] newChangePoints, double precision){
        int currentLength=currentChangePoints.length;
        int addLength=newChangePoints.length;
        int newLength = currentLength+addLength;

//       Takes the currentChangePoints and the newChangePoints into getData() but ordered and the order in getOrder()
//        assumes that order numbers 0..currentLength are the positions of currentChangePoints
//        and currentLength..newLength  - currentLentgth are the positions of newChangePoints

        QuadupleGP tempQuad= sortUpdate(currentChangePoints, newChangePoints);
//          index=tempPair.getOrder();
//        totalChangePoints=tempPair.getData();
//      This is a silly thing to only compute the Q.matrix for the neighbors of the newChangePoints
//      Takes the positions of the newData in the new complete sorted data and adds the positions of neighbors


//      Takes the positions of the newData in the new complete sorted data and adds the positions of neighbors
        int [] NeighborsIndex =Neighbors(tempQuad.getPositionNew(), newLength);

//      Retrieves the positions indicated in NeigborsIndex from the complete sorted data
        no.uib.cipr.matrix.DenseVector tempData = new DenseVector(SubsetData(tempQuad.getData(),NeighborsIndex));

        SymmTridiagMatrix Q = getQmatrix(precision, tempData);



//        Retrieves the positions indicated in NeighborsIndex from the getOrder
        int [] NeighborsOriginal= SubsetData(tempQuad.getOrder(),NeighborsIndex);
//      Generates two arrays: one with the positions of newData and other with the positions of OldData in the Neighbors data TempData
        PairIndex Indicators = SubIndex(NeighborsOriginal,currentLength,addLength);

        UpperSPDBandMatrix varf = new UpperSPDBandMatrix(Matrices.getSubMatrix(Q, Indicators.getOrderNew(), Indicators.getOrderNew()),1);
        BandCholesky U1 = new BandCholesky(addLength,1,true);

        U1.factor(varf);



        DenseVector part = new DenseVector(addLength);
        int [] GPpositions =SubsetData(NeighborsOriginal,Indicators.getOrderOld());
        DenseVector currentGPvalues= new DenseVector(currentGPvalues1);
        DenseVector currentGPneighbors = new DenseVector(SubsetData(currentGPvalues,GPpositions));


        Matrix first = Matrices.getSubMatrix(Q,Indicators.getOrderNew(),Indicators.getOrderOld());

        first.mult(-1,currentGPneighbors,part);

        DenseVector mean = new DenseVector(getMultiNormalMean(part, U1.getU()));

        double [] addGPvalues = getMultiNormal(mean,U1.getU()).getData();

        return addGPvalues;

    }

//    public PairGP sortNew(double [] newPoints){
//        double [] newPoints2= new double[newPoints.length];
//        int [] newOrder=new int[newPoints.length];

//
//        double pivot=newPoints[0];
//        for (int j=0; j<newPoints.length;j++){
//
//        }
//                                                          newPoints2;
//
//
//        return new PairGP(newPoints2,newOrder);
//    }

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


        int [] NeighborsIndex =Neighbors(tempQuad.getPositionNew());


//      Retrieves the positions indicated in NeigborsIndex from the complete sorted data
        DenseVector tempData = new DenseVector(SubsetData(tempQuad.getData(),NeighborsIndex));

        SymmTridiagMatrix Q = getQmatrix(precision, tempData);


//        Retrieves the positions indicated in NeighborsIndex from the getOrder
//        int [] NeighborsOriginal= SubsetData(tempQuad.getOrder(),NeighborsIndex);

//        System.err.println("neighborsOriginal"+NeighborsOriginal[0]+"b:"+NeighborsOriginal[1]+"c"+NeighborsOriginal[2]);
//

        double part=0.0;
        double varf=1.0;
        if (NeighborsIndex.length==3){
           varf = Q.get(1,1);


            part = -currentGPvalues.get(NeighborsIndex[0])*Q.get(1,0)-currentGPvalues.get(NeighborsIndex[1])*Q.get(1,2);

        }
        if (NeighborsIndex.length==2){
           varf = Q.get(0,0);
           part = -currentGPvalues.get(NeighborsIndex[0])*Q.get(0,1);
        }
        double res =(part/varf)+(MathUtils.nextGaussian()/Math.sqrt(varf));

        return new Trip1GP(res,tempQuad.getOrder(),tempQuad.getPositionNew());
        }

    public double[] getGPvalue(double [] currentChangePoints, DenseVector currentGPvalues, double newChangePoints, double precision){

        int currentLength=currentChangePoints.length;
        int addLength=1;
//        int newLength = currentLength+addLength;
        Quaduple1GP tempQuad= sortUpdate(currentChangePoints, newChangePoints);
        int [] NeighborsIndex =Neighbors(tempQuad.getPositionNew());
        DenseVector tempData = new DenseVector(SubsetData(tempQuad.getData(),NeighborsIndex));


        SymmTridiagMatrix Q = getQmatrix(precision, tempData);

        double part=0.0;
        double varf=1.0;
        if (NeighborsIndex.length==3){
            varf = Q.get(1,1);

            part = -currentGPvalues.get(NeighborsIndex[0])*Q.get(1,0)-currentGPvalues.get(NeighborsIndex[1])*Q.get(1,2);

        }
        if (NeighborsIndex.length==2){
            varf = Q.get(0,0);
            part = -currentGPvalues.get(NeighborsIndex[0])*Q.get(0,1);
        }
        double res =(part/varf)+(MathUtils.nextGaussian()/Math.sqrt(varf));
        double logpr=0.5*Math.log(varf)-0.5*varf*Math.pow(res-(part/varf),2);
        double [] end=new double[3];
        end[0]=res;
        end[1]=logpr;
        end[2]=tempQuad.getPositionNew();
//        System.err.println("here"+end[2]);
        return end;
    }

    public double[] getGPvalueroot(double [] currentChangePoints, DenseVector currentGPvalues, double newChangePoints, double precision){
//Here, we know that the newChangePoint is larger than all the currentChangePoints, so we only need one currentChangePoint
//       and one currentGPvalues (the last one)

//        int currentLength=currentChangePoints.length;
//        int addLength=1;
//        int newLength = currentLength+addLength;
//        Quaduple1GP tempQuad= sortUpdate(currentChangePoints, newChangePoints);
//        int [] NeighborsIndex =Neighbors(tempQuad.getPositionNew());
//        DenseVector tempData = new DenseVector(SubsetData(tempQuad.getData(),NeighborsIndex));
        double [] tempD=new double[2];
        double part=0.0;
        double varf=1.0;
        if (currentChangePoints[currentChangePoints.length-1]<newChangePoints){
        tempD[0]=currentChangePoints[currentChangePoints.length-1];
        tempD[1]=newChangePoints;
            DenseVector tempData=new DenseVector(tempD);
//        System.err.println(tempData);
            SymmTridiagMatrix Q = getQmatrix(precision, tempData);


            varf = Q.get(1,1);
            part = -currentGPvalues.get(currentGPvalues.size()-1)*Q.get(0,1);
        }
        else {
            tempD[1]=currentChangePoints[currentChangePoints.length-1];
            tempD[0]=newChangePoints;
            DenseVector tempData=new DenseVector(tempD);
//        System.err.println(tempData);
            SymmTridiagMatrix Q = getQmatrix(precision, tempData);


            varf = Q.get(0,0);
            part = -currentGPvalues.get(currentGPvalues.size()-1)*Q.get(0,1);


        }

        double res =(part/varf)+(MathUtils.nextGaussian()/Math.sqrt(varf));
        double logpr=0.5*Math.log(varf)-0.5*varf*Math.pow(res-(part/varf),2);
        double [] end=new double[2];
        end[0]=res;
        end[1]=logpr;
        return end;
    }

    public double getDensity(double [] currentChangePoints, DenseVector currentGPvalues, double newChangePoints, double precision,int who){
        double [] temp;
        SymmTridiagMatrix Q;
        double part=0.0;
        double varf=1.0;
//                           System.err.println("who"+who);

        if (who==0){
            temp=new double[2];
            temp[0]=newChangePoints;
            temp[1]=currentChangePoints[1];
            DenseVector tempData = new DenseVector(temp);


            Q = getQmatrix(precision, tempData);
            varf = Q.get(0,0);
            part = -currentGPvalues.get(1)*Q.get(0,1);
            }

        if (who==currentChangePoints.length-1){
            temp=new double[2];
            temp[0]=currentChangePoints[who-1];
            temp[1]=newChangePoints;
            DenseVector tempData = new DenseVector(temp);
            Q = getQmatrix(precision, tempData);
            varf = Q.get(1,1);
            part = -currentGPvalues.get(0)*Q.get(0,1);
        }
        if (who>0 & who<currentChangePoints.length-1){
//            System.err.println("whooo");
            temp=new double[3];
            temp[0]=currentChangePoints[who-1];
            temp[1]=newChangePoints;
            temp[2]=currentChangePoints[who+1];
            DenseVector tempData = new DenseVector(temp);
            Q = getQmatrix(precision, tempData);
            varf = Q.get(1,1);
             part = -currentGPvalues.get(who-1)*Q.get(1,0)-currentGPvalues.get(who+1)*Q.get(1,2);


        }
//        Quaduple1GP tempQuad= sortUpdate(currentChangePoints, newChangePoints);
//        int [] NeighborsIndex =Neighbors(tempQuad.getPositionNew());
        double res =currentGPvalues.get(who);
        double logpr=0.5*Math.log(varf)-0.5*varf*Math.pow(res-(part/varf),2);
        return logpr;
    }

    public double sigmoidal (double x){
       return 1/(1+Math.exp(-x));
   }

//
    public void addOnePoint(double gpval, double location, int positions, int where){
       GPcounts.setParameterValue(where,GPcounts.getParameterValue(where)+1.0);
        popSizeParameter.addDimension(positions,gpval);
       changePoints.addDimension(positions, location);
       GPtype.addDimension(positions,-1.0);
       coalfactor.addDimension(positions,GPcoalfactor[where]);
    }

     public void delOnePoint(int positions, int locIntervals){
       popSizeParameter.removeDimension(positions);
       changePoints.removeDimension(positions);
       GPtype.removeDimension(positions);
       coalfactor.removeDimension(positions);
       GPcounts.setParameterValue(locIntervals,GPcounts.getParameterValue(locIntervals)-1.0);
    }

    public int wherePoint(double point, int x, double lower) {
        int z=0;
        int found=0;
        if (x==-1){
            z=0;
        }
        if (x>=0){
            z=CoalPosIndicator[x]+1;
        }
            while (found==0){
                if (point<=lower+GPvalue.getInterval(z)){
                    found=1;
                }  else {
                    lower+=GPvalue.getInterval(z);
                    z++;
                }
            }

        return z;
    }

    public int searchPos(double[] dataFind, double oldValue, double newValue, int currentPosition){
        boolean found=false;
        if (oldValue<newValue){
//            look back
            while (found==false){
                currentPosition--;
                if (dataFind[currentPosition]<=oldValue){
                    found=true;
                }
            }
        }  else {
            while (found==false){
                currentPosition++;
                if(dataFind[currentPosition]>=oldValue){
                    found=true;
                }
            }
        }
        return currentPosition;
    }

    public void writeChain(double [] current,String fileName){
        try {
            BufferedWriter writer = new BufferedWriter(new FileWriter(fileName,true));
            for (int i=0;i<current.length;i++){
                writer.write(current[i]+" ");
            }
            writer.newLine();
            writer.close();
        }
        catch(java.io.IOException ioe){
                System.err.println("IOException:"+ ioe.getMessage());
            }
        }


       public void numberThinned(double [] currentChangePoints, DenseVector currentPopSize, double currentPrecision) {
           numberPoints=CoalCounts.getSize();
           double [] addPoints=new double[numberPoints];
           int [] indAPoint = new int[numberPoints];
           int [] indicator =new int[numberPoints];
           double lowL = 0.0;
           double uppL =0.0;
           TripGP tempG;
           double accept=0.0;
           double ltime=0.0;
           int where=0;
           int who=0;
           int k=0;
           int k1=0;
           int fix=0;


//         Proposes to add/delete and proposes location uniformly in each intercoalescent interval
       for (int j=0; j<numberPoints; j++){
           uppL+=GPvalue.getGPCoalInterval(j);


           if (MathUtils.nextDouble()<0.5) {
               indicator[j]=1;
               addPoints[k]=MathUtils.uniform(lowL,uppL);
               indAPoint[k]=wherePoint(addPoints[k],j-1,lowL);
               k++;
           }
            lowL=uppL;
       }

           double [] addPoints2 = new double[k];

           System.arraycopy(addPoints,0,addPoints2,0,k);
           tempG = getGPvalues(currentChangePoints, currentPopSize, addPoints2,currentPrecision);
           for (int i=0; i<numberPoints;i++){
               if (indicator[i]==1){
                   accept=lambdaBoundParameter.getParameterValue(0)*GPvalue.getGPCoalInterval(i)*GPcoalfactor[indAPoint[k1]]*sigmoidal(-tempG.getData()[k1])/(CoalCounts.getParameterValue(i)+1);
               if (MathUtils.nextDouble()<accept) {
                   addOnePoint(tempG.getData()[k1], addPoints2[k1], tempG.getNewOrder()[k1] + fix, indAPoint[k1]);
                   CoalCounts.setParameterValue(i,CoalCounts.getParameterValue(i)+1);
               }  else {
                   fix--;}
                   k1++;

               }   else {
                   double temp=CoalCounts.getParameterValue(i);
                   if (temp>1){
//                         it's 0 inclusive, no need to -1
                       int rdm = MathUtils.nextInt((int) temp-1);
                       who=where+rdm;
                       accept=temp/(lambdaBoundParameter.getParameterValue(0)*GPvalue.getGPCoalInterval(i)*sigmoidal(-popSizeParameter.getParameterValue(who))*coalfactor.getParameterValue(who));
                            if (MathUtils.nextDouble()<accept){
                                if (i>0) {ltime=CoalTime[i-1];} else {ltime=0.0;}
                                int locationIntervals=wherePoint(changePoints.getParameterValue(who),i-1, ltime);
                                delOnePoint(who+1,locationIntervals);
                                CoalCounts.setParameterValue(i,CoalCounts.getParameterValue(i)-1);
                                fix--;

                }
               }
                 if (temp==1){
                     if (i>0) {ltime=CoalTime[i-1];} else {ltime=0.0;}
                     int locationIntervals=wherePoint(changePoints.getParameterValue(where), i - 1,  ltime);
                     delOnePoint(where+1,locationIntervals);
                     CoalCounts.setParameterValue(i,0);
                     fix--;
                 }
               }
               where+=CoalCounts.getParameterValue(i)+1;
           }

       }


    public void locationThinned(double [] currentChangePoints, DenseVector currentPopSize, double currentPrecision) {
//    Sample the interval
        double where;
        double lowL = 0.0;
        double [] toAdd;
        double [] toDel;
        double uppL =0.0;
        double cum=0.0;
        boolean got=false;
        double Tlen=0.0;
        double numberPoints=CoalCounts.getSize();
        for (int j=0; j<numberPoints; j++){
            if (CoalCounts.getParameterValue(j)>0) {
          Tlen+=GPvalue.getGPCoalInterval(j);
            }
        }
//          Need to only consider the intevals with latent points
        int j=0;
        int [] ind;
        int [] ind2;
        double [] toDelG;
        int position =0;
        int oldPos=0;
        double accept;
        if (Tlen>0){
        where= MathUtils.uniform(0,Tlen);
        while (!got) {
            cum+=GPvalue.getGPCoalInterval(j);
            if (CoalCounts.getParameterValue(j)>0){
                 uppL+=GPvalue.getGPCoalInterval(j);
                 if (where>lowL & where<=uppL){
                     got=true;
                 }

                lowL=uppL;
            }
        position+=CoalCounts.getParameterValue(j)+1;
            j++;
        }
        j--;
        position-=CoalCounts.getParameterValue(j)+1;
        Trip1GP tempG;
        int m=(int) CoalCounts.getParameterValue(j);
        toAdd=new double[m];
        toDelG=new double[m];
        ind= new int[m];
        toDel= new double[m];
        ind2= new int[m];

        for (int i=0; i<m; i++){
            toAdd[i]=MathUtils.uniform(cum-GPvalue.getGPCoalInterval(j),cum);
//            System.err.println("add"+toAdd[i]+"from:"+(cum-GPvalue.getGPCoalInterval(j))+"check"+CoalTime[j]);
            ind[i]=wherePoint(toAdd[i],j-1,cum-GPvalue.getGPCoalInterval(j));
            toDel[i]=currentChangePoints[position+i];
            toDelG[i]=currentPopSize.get(position+i);
//            System.err.println("verify:"+toDel[i]+"be replaced by"+toAdd[i]+"in?"+ind[i]+"interval:"+(cum-GPvalue.getGPCoalInterval(j))+"and"+cum);
            ind2[i]=wherePoint(changePoints.getParameterValue(position+i),j-1,cum-GPvalue.getGPCoalInterval(j));
//            System.err.println(i+": to add:"+toAdd[i]+"and pos"+position);
        }

//


//        Arrays.sort(index,toAdd);
        for (int i=0; i<m; i++){
            tempG = getGPvalues(currentChangePoints, currentPopSize, toAdd[i],currentPrecision);


//            tempG=getGPvalues(currentChangePoints,currentPopSize,toAdd[i],currentPrecision);
//            System.err.println(tempG.getNewOrder()+" as new order");
//            System.err.println("gp fac old:"+GPcoalfactor[ind2[i]]+" gp factor new"+GPcoalfactor[ind[i]]);
//            System.err.println("gp val old:"+toDelG[i]+" gp val new:"+tempG.getData());
//
            accept=(GPcoalfactor[ind[i]]*sigmoidal(-tempG.getData()))/(GPcoalfactor[ind2[i]]*(sigmoidal(-toDelG[i])));

            if (MathUtils.nextDouble()<accept){
//                      System.err.println("here2");
                addOnePoint(tempG.getData(),toAdd[i],tempG.getNewOrder(),ind[i]);
                currentChangePoints=(double []) changePoints.getParameterValues();
                currentPopSize=new DenseVector(popSizeParameter.getParameterValues());

                oldPos=searchPos(currentChangePoints,toDel[i],toAdd[i],tempG.getNewOrder());
//
//                System.err.println("removes in location"+toDel[i]+"pos"+oldPos+"value?"+changePoints.getParameterValue(oldPos)+"from"+ind2[i]);
                delOnePoint(oldPos+1,ind2[i]);
                currentChangePoints = (double []) changePoints.getParameterValues();
                currentPopSize=new DenseVector(popSizeParameter.getParameterValues());
            }
        }
    }
    }

    public double forLikelihood(double [] Gvalues, Parameter Gtype){
        double loglik=0.0;
        for (int j=0;j<Gvalues.length;j++){

          loglik-=Math.log(1+Math.exp(-Gtype.getParameterValue(j)*Gvalues[j]));
          }
        return loglik;
    }

    public void sliceSampling(double [] currentChangePoints, DenseVector currentPopSize, double currentPrecision){
        double theta;
        double thetaPrime=0.0;
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
//        System.err.println("loglik"+loglik);
        double loglik2=0.0;

        while (keep==1){
            thetaPrime=MathUtils.uniform(thetaMin,thetaMax);
//            System.err.println("theta Prime"+thetaPrime+"theta"+theta+"thetaMin"+thetaMin+"thetaMax"+thetaMax);
            proposal.zero();
            proposal.add(Math.sin(thetaPrime),v1);
            proposal.add(Math.cos(thetaPrime),v2);
            double [] popSize2 = proposal.getData();
            loglik2=forLikelihood(popSize2,GPtype);
//            System.err.println("loglik2"+loglik2);

            if (loglik2>loglik) {keep=2;}
            else {
                if (thetaPrime<theta) {thetaMin=thetaPrime;} else {thetaMax=thetaPrime;}
            }
        }
        for (int j=0; j<popSizeParameter.getSize();j++)
            popSizeParameter.setParameterValue(j, proposal.get(j));
//            popSizeParameter.setParameterValueQuietly(j,proposal.get(j));

//        ((Parameter.Abstract) popSizeParameter).fireParameterChangedEvent();

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
    public double doOperation() {

        System.err.println("Runs my operator");

        double currentPrecision = this.precisionParameter.getParameterValue(0);
        DenseVector currentPopSize = new DenseVector(popSizeParameter.getParameterValues());
        double [] currentChangePoints = this.changePoints.getParameterValues();

        currentQ=getQmatrix(currentPrecision,currentChangePoints);
        double currentQuadratic = getQuadraticForm(currentQ, currentPopSize);
        double newprecision=getNewPrecision(currentPrecision,currentQuadratic);
//         Gibbs sample new precision
        precisionParameter.setParameterValue(0,newprecision);
        currentPrecision=this.precisionParameter.getParameterValue(0);


        double currentLambda = this.lambdaBoundParameter.getParameterValue(0);
        getNewUpperBound(currentLambda);
        currentLambda=this.lambdaBoundParameter.getParameterValue(0);
//        System.err.println("before thinned:"+changePoints.getParameterValue(75));
//
        numberThinned(currentChangePoints, currentPopSize, currentPrecision);
//        System.err.println("after thinned:"+changePoints.getParameterValue(75));
//
        DenseVector currentPopSize1 = new DenseVector(popSizeParameter.getParameterValues());
        double [] currentChangePoints1 = this.changePoints.getParameterValues();
//        System.err.println("finished number"+popSizeParameter.getSize());

        locationThinned(currentChangePoints1,currentPopSize1,currentPrecision);

//        System.err.println("finished location"+popSizeParameter.getSize());

        DenseVector currentPopSize2 = new DenseVector(popSizeParameter.getParameterValues());
//        System.err.println("2:"+currentPopSize2);

        double [] currentChangePoints2 = this.changePoints.getParameterValues();
//        System.exit(-1);

        sliceSampling(currentChangePoints2, currentPopSize2, currentPrecision);
        numPoints.setParameterValue(0,popSizeParameter.getSize());

        double [] currentPopSize3 = this.popSizeParameter.getParameterValues();

        int sumcoal=0;
        int sumlat=0;

        for (int j=0;j<changePoints.getSize();j++){
            if( GPtype.getParameterValue(j)==1) {sumcoal++;}

        }
        for (int j=0;j<CoalCounts.getSize();j++){
           sumlat+=CoalCounts.getParameterValue(j);
        }
        if (sumcoal!=CoalCounts.getSize()){ System.err.println("WARNING CONSISTENCY 1");   }
        if (sumlat!=(changePoints.getSize()-CoalCounts.getSize())){ System.err.println("WARNING CONSISTENCY 2"+sumlat+"and changePts size"+changePoints.getSize()); }

//        System.err.println("check values 72:"+changePoints.getParameterValue(72)+"type"+GPtype.getParameterValue(72));
//        System.err.println("check values 73:"+changePoints.getParameterValue(73)+"type"+GPtype.getParameterValue(73));

//        System.err.println("check values 74:"+changePoints.getParameterValue(74)+"type"+GPtype.getParameterValue(74)+"counts"+GPcounts.getParameterValue(114));
//        System.err.println("check values 75:"+changePoints.getParameterValue(75)+"type"+GPtype.getParameterValue(75)+"counts"+GPcounts.getParameterValue(115));
//        System.err.println("check values 76:"+changePoints.getParameterValue(76)+"type"+GPtype.getParameterValue(76)+"counts"+GPcounts.getParameterValue(116));
//        System.err.println("check values 3:"+changePoints.getParameterValue(3)+"type"+GPtype.getParameterValue(3)+"counts"+GPcounts.getParameterValue(16));
//        System.err.println("check values 4:"+changePoints.getParameterValue(4)+"type"+GPtype.getParameterValue(4)+"counts"+GPcounts.getParameterValue(18));
//        System.err.println("check values 5:"+changePoints.getParameterValue(5)+"type"+GPtype.getParameterValue(5)+"counts"+GPcounts.getParameterValue(19));
//
//        writeChain(currentChangePoints2,"locations.txt");
//        writeChain(currentPopSize3,"gvalues.txt");



//
////
//        System.err.println("type after"+GPtype.getSize());

//
//
////        ArrayList<ComparableDouble> times = new ArrayList<ComparableDouble>();
////                ArrayList<Integer> childs = new ArrayList<Integer>();
////                collectAllTimes(tree, root, exclude, times, childs);
////                int[] indices = new int[times.size()];
//
//
//
////        double currentPrecision = precisionParameter.getParameterValue(0);
//        double proposedPrecision = this.getNewPrecision(currentPrecision, scaleFactor);
//
////        double currentLambda = this.lambdaParameter.getParameterValue(0);
//        double proposedLambda = currentLambda;
//
//        precisionParameter.setParameterValue(0, proposedPrecision);
//        lambdaParameter.setParameterValue(0, proposedLambda);
//
//        DenseVector currentGamma = new DenseVector(GPvalue.getPopSizeParameter().getParameterValues());
//        DenseVector proposedGamma;
//
////        SymmTridiagMatrix currentQ = GPvalue.getStoredScaledWeightMatrix(currentPrecision, currentLambda);
////        SymmTridiagMatrix proposedQ = GPvalue.getScaledWeightMatrix(proposedPrecision, proposedLambda);
//
//
//
//
////        double[] wNative = gmrfField.getSufficientStatistics();
//
////        UpperSPDBandMatrix forwardQW = new UpperSPDBandMatrix(proposedQ, 1);
////        UpperSPDBandMatrix backwardQW = new UpperSPDBandMatrix(currentQ, 1);
////
////        BandCholesky forwardCholesky = new BandCholesky(wNative.length, 1, true);
////        BandCholesky backwardCholesky = new BandCholesky(wNative.length, 1, true);
////
////        DenseVector diagonal1 = new DenseVector(fieldLength);
////        DenseVector diagonal2 = new DenseVector(fieldLength);
////        DenseVector diagonal3 = new DenseVector(fieldLength);
////
////        DenseVector modeForward = newtonRaphson(wNative, currentGamma, proposedQ.copy());
////
////        for (int i = 0; i < fieldLength; i++) {
////            diagonal1.set(i, wNative[i] * Math.exp(-modeForward.get(i)));
////            diagonal2.set(i, modeForward.get(i) + 1);
////
////            forwardQW.set(i, i, diagonal1.get(i) + forwardQW.get(i, i));
////            diagonal1.set(i, diagonal1.get(i) * diagonal2.get(i) - 1);
////        }
////
////        forwardCholesky.factor(forwardQW.copy());
////
////        DenseVector forwardMean = getMultiNormalMean(diagonal1, forwardCholesky);
////
////        DenseVector stand_norm = new DenseVector(zeros);
////
////        for (int i = 0; i < zeros.length; i++)
////            stand_norm.set(i, MathUtils.nextGaussian());
////
////        proposedGamma = getMultiNormal(stand_norm, forwardMean, forwardCholesky);
////
////
////        for (int i = 0; i < fieldLength; i++)
////            popSizeParameter.setParameterValueQuietly(i, proposedGamma.get(i));
////
////        ((Parameter.Abstract) popSizeParameter).fireParameterChangedEvent();
////
//
//        double hRatio = 0;
//
////        diagonal1.zero();
////        diagonal2.zero();
////        diagonal3.zero();
////
////        DenseVector modeBackward = newtonRaphson(wNative, proposedGamma, currentQ.copy());
////
////        for (int i = 0; i < fieldLength; i++) {
////            diagonal1.set(i, wNative[i] * Math.exp(-modeBackward.get(i)));
////            diagonal2.set(i, modeBackward.get(i) + 1);
////
////            backwardQW.set(i, i, diagonal1.get(i) + backwardQW.get(i, i));
////            diagonal1.set(i, diagonal1.get(i) * diagonal2.get(i) - 1);
////        }
////
////        backwardCholesky.factor(backwardQW.copy());
////
////        DenseVector backwardMean = getMultiNormalMean(diagonal1, backwardCholesky);
////
////        for (int i = 0; i < fieldLength; i++) {
////            diagonal1.set(i, currentGamma.get(i) - backwardMean.get(i));
////        }
////
////        backwardQW.mult(diagonal1, diagonal3);
////
////        // Removed 0.5 * 2
////        hRatio += logGeneralizedDeterminant(backwardCholesky.getU()) - 0.5 * diagonal1.dot(diagonal3);
////        hRatio -= logGeneralizedDeterminant(forwardCholesky.getU() ) - 0.5 * stand_norm.dot(stand_norm);
//
//
////        return hRatio;
////        System.err.println("Prueba");
////        System.exit(-1);
        return 0.0;
    }

    //MCMCOperator INTERFACE
  // This is the only part where GPSBUOperateroParser is used
    public final String getOperatorName() {
        return GaussianProcessSkytrackBlockUpdateOperatorParser.BLOCK_UPDATE_OPERATOR;
    }


    public  double getTemperature() {
        return 0.0;
    }


    public double getCoercableParameter() {
//        return Math.log(scaleFactor);
        return Math.sqrt(scaleFactor - 1);
    }

    public int getStepCount() {
        return 1;
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
