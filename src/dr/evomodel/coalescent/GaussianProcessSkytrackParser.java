package dr.evomodel.coalescent;

import dr.math.MathUtils;
import dr.stats.DiscreteStatistics;
import no.uib.cipr.matrix.*;

import java.io.BufferedReader;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.Arrays;

/**
 * Created by mkarcher on 6/22/16.
 */
public class GaussianProcessSkytrackParser {
    private static final int changepointsIndex = 5;
    private static final int GvaluesIndex = 6;
    private static final int lambdaBoundIndex = 7;
    private static final int precisionIndex = 8;
    private static final int tmrcaIndex = 10;

    public static void main(String[] args) {
//        int [] indices = {0, 5, 9};
//        int total = 10;
//        int [] neighbs = neighbors(indices, total);
//        System.out.println(Arrays.toString(neighbs));

//        double [] oldData = {1.0, 2.0, 3.0};
//        double [] newData = {1.5, 3.5};
//
//        QuadupleGP result = sortUpdate(oldData, newData);
//
//        System.out.println("Data: " + Arrays.toString(result.getData()));
//        System.out.println("Order: " + Arrays.toString(result.getOrder()));
//        System.out.println("New Positions: " + Arrays.toString(result.getPositionNew()));
//        System.out.println("Old Positions: " + Arrays.toString(result.getPositionOld()));

        int NGRID = 101;
        double quantile = 0.5;
        String filePathName = "examples/hcvNew2small.log";
        CSVstats stats = parseCSV(filePathName, 3);
        //System.out.println(Arrays.toString(stats.precisions));
        double tmrca = DiscreteStatistics.quantile(quantile, stats.tmrcas);
        System.out.println(tmrca);

        double[] grid = new double[NGRID];
        double step = tmrca / (NGRID-1);
        grid[0] = 0.001;
        for (int i = 1; i < grid.length; i++) {
            grid[i] = grid[i-1] + step;
        }

        double[] result = gpPosterior(stats, grid, stats.getSize() - 1);

        double[] result2 = new double[result.length];
        for (int i = 0; i < result.length; i++) {
            result2[i] = 1/result[i];
        }

        System.out.println(Arrays.toString(grid));
        System.out.println(Arrays.toString(result2));
    }

    public static CSVstats parseCSV(String filename, int skip) {
        CSVstats result = null;

        try {
            BufferedReader br = new BufferedReader(new FileReader(filename));

            String line;
            String[] fieldsArray;
            int lineNum = 0;

            ArrayList<ArrayList<Double>> changepoints = new ArrayList<ArrayList<Double>>();
            ArrayList<ArrayList<Double>> Gvalues = new ArrayList<ArrayList<Double>>();
            ArrayList<Double> lambdas = new ArrayList<Double>();
            ArrayList<Double> precisions = new ArrayList<Double>();
            ArrayList<Double> tmrcas = new ArrayList<Double>();

            while ((line = br.readLine()) != null) {
                if (lineNum >= skip) {
                    fieldsArray = line.split("\t");

                    changepoints.add(parseListStr(fieldsArray[changepointsIndex]));
                    Gvalues.add(parseListStr(fieldsArray[GvaluesIndex]));
                    lambdas.add(new Double(fieldsArray[lambdaBoundIndex]));
                    precisions.add(new Double(fieldsArray[precisionIndex]));
                    tmrcas.add(new Double(fieldsArray[tmrcaIndex]));
                }
                lineNum += 1;
            }

            result = new CSVstats(changepoints, Gvalues, lambdas, precisions, tmrcas);
        } catch (Exception e) {
            e.printStackTrace();
        }

        return result;
    }

    public static double[] gpPosterior(CSVstats stats, double[] grid, int iter) {
        double[] changePoints = stats.changepoints.get(iter);
        DenseVector values = new DenseVector(stats.Gvalues.get(iter));
        double precision = stats.precisions[iter];

        TripGP res1 = getGPvalues(changePoints, values, grid, precision);

        double[] result1 = sigmoidal(res1.getData());

        double[] result = new double[result1.length];
        for (int i = 0; i < result1.length; i++) {
            result[i] = result1[i] * stats.lambdas[iter];
        }

        return result;
    }

    public static TripGP getGPvalues(double[] currentChangePoints, DenseVector currentGPvalues,
                                     double[] newChangePoints, double precision) {
        int currentLength=currentChangePoints.length;
        int addLength=newChangePoints.length;
        int newLength = currentLength+addLength;

//       Takes the currentChangePoints and the newChangePoints into getData() but ordered and the order in getOrder()
//        assumes that order numbers 0..currentLength are the positions of currentChangePoints
//        and currentLength..newLength  - currentLentgth are the positions of newChangePoints

        QuadupleGP tempQuad = sortUpdate(currentChangePoints, newChangePoints);
//          index=tempPair.getOrder();
//        totalChangePoints=tempPair.getData();
//      This is a silly thing to only compute the Q.matrix for the neighbors of the newChangePoints
//      Takes the positions of the newData in the new complete sorted data and adds the positions of neighbors



//      Takes the positions of the newData in the new complete sorted data and adds the positions of neighbors
        int [] NeighborsIndex = neighbors(tempQuad.getPositionNew(), newLength);

//      Retrieves the positions indicated in NeigborsIndex from the complete sorted data
        DenseVector tempData = new DenseVector(SubsetData(tempQuad.getData(),NeighborsIndex));

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
        DenseVector currentGPneighbors = new DenseVector(SubsetData(currentGPvalues,GPpositions));

        Matrix first = Matrices.getSubMatrix(Q,Indicators.getOrderNew(),Indicators.getOrderOld());

        first.mult(-1,currentGPneighbors,part);

        DenseVector mean = new DenseVector(getMultiNormalMean(part, U1.getU()));

        double [] addGPvalues = getMultiNormal(mean,U1.getU()).getData();

        return new TripGP(addGPvalues,tempQuad.getOrder(),tempQuad.getPositionNew());
    }

    private static DenseVector getMultiNormalMean(DenseVector CanonVector, UpperTriangBandMatrix CholeskyUpper) {
        DenseVector tempValue = new DenseVector(CanonVector.size());
        DenseVector Mean = new DenseVector(CanonVector.size());

        // UpperTriangBandMatrix CholeskyUpper = Cholesky.getU();

        // Assume Cholesky factorization of the precision matrix Q = LL^T
        // 1. Solve L\omega = b
        // System.err.println(CanonVector.size()+"and"+tempValue.size());
        CholeskyUpper.transSolve(CanonVector, tempValue);

        // 2. Solve L^T \mu = \omega
        CholeskyUpper.solve(tempValue, Mean);

        return Mean;
    }

    private static DenseVector getMultiNormal(DenseVector Mean, UpperTriangBandMatrix CholeskyUpper) {
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

    //  Returns an array with the positions indicated in Index
    private static double [] SubsetData(double [] Data, int [] Index) {
        double [] Res= new double [Index.length];
        for (int j=0;j<Index.length;j++){
            Res[j]=Data[Index[j]];
        }
        return Res;
    }

    private static double [] SubsetData(DenseVector Data, int [] Index) {
        double [] Res= new double [Index.length];
        for (int j=0;j<Index.length;j++){
            Res[j]= Data.get(Index[j]);
        }
        return Res;
    }
    private static int [] SubsetData(int [] Data, int [] Index) {
        int [] Res= new int [Index.length];
        for (int j=0;j<Index.length;j++){
            Res[j]=Data[Index[j]];
        }
        return Res;
    }

    private static ArrayList<Double> parseListStr(String arg) {
        ArrayList<Double> vals = new ArrayList<Double>();

        String[] valStrs = arg.replaceAll("\\{|}", "").split(",");
        for (String valStr : valStrs) {
            vals.add(new Double(valStr));
        }

        return vals;
    }

    private static double sigmoidal(double x) {
        return 1/(1+Math.exp(-x));
    }

    private static double[] sigmoidal(double[] xs) {
        double[] result = new double[xs.length];
        for (int i = 0; i < result.length; i++) {
            result[i] = sigmoidal(xs[i]);
        }
        return result;
    }

    //Returns the index (position) of newData + Neighbors in the Ordered List
    /*private static int [] neighbors(int [] indexNew, int numberTotalData){
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
    }*/

    private static int [] neighbors(int [] indexNew, int numberTotalData) {
        int [] ns = new int[numberTotalData];
        int index, k = 0;
        int [] result1 = new int[numberTotalData];
        int [] result;

        for (int i = 0; i < indexNew.length; i++) {
            index = indexNew[i];

            if (index-1 > 0) {
                ns[index-1]++;
            }

            ns[index]++;

            if (index+1 < numberTotalData) {
                ns[index+1]++;
            }
        }

        for (int i = 0; i < numberTotalData; i++) {
            if (ns[i] > 0) {
                result1[k] = i;
                k++;
            }
        }

        result = new int[k];
        System.arraycopy(result1, 0, result, 0, k);

        return result;
    }

    // Assumes the input vector x is ordered
    private static SymmTridiagMatrix getQmatrix(double precision, DenseVector x) {
        SymmTridiagMatrix res;
        double trick = 0;
        double[] offdiag = new double[x.size() - 1];
        double[] diag = new double[x.size()];

        for (int i = 0; i < x.size() - 1; i++) {
            offdiag[i] = precision*(-1.0 / (x.get(i+1)-x.get(i)));
            if (i< x.size()-2){
                diag[i+1]= -offdiag[i]+precision*(1.0/(x.get(i+2)-x.get(i+1))+trick);
            }
        }
        // Diffuse prior correction - intrinsic
        // Take care of the endpoints
        diag[0] = -offdiag[0]+precision*trick;

        diag[x.size() - 1] = -offdiag[x.size() - 2]+precision*(trick);
        res = new SymmTridiagMatrix(diag, offdiag);
        return res;
    }

    private static SymmTridiagMatrix getQmatrix(double precision, double[] x) {
        SymmTridiagMatrix res;
        double trick = 0.00000000001;
        double[] offdiag = new double[x.length - 1];
        double[] diag = new double[x.length];

        for (int i = 0; i < x.length - 1; i++) {
            offdiag[i] = precision*(-1.0 / (x[i+1]-x[i]));

            if (i< x.length-2){
                diag[i+1]= -offdiag[i]+precision*(1.0/(x[i+2]-x[i+1])+trick);
            }
        }
        // Diffuse prior correction - intrinsic
        // Take care of the endpoints
        diag[0] = -offdiag[0]+precision*trick;

        diag[x.length - 1] = -offdiag[x.length - 2]+precision*(trick);
        res = new SymmTridiagMatrix(diag, offdiag);
        return res;
    }

    private static QuadupleGP sortUpdate(double [] sortedData, double [] newData){
        // note that sortedData and newData are already ordered (minimum to maximum)
        // and last(sortedData) > last(newData)
        int newLength = sortedData.length + newData.length;
        double [] newList = new double [newLength];
        int [] newOrder = new int [newLength];

        // indexNew contains the index where the newData is stored (index ordered) in newList
        // indexOld contains the index where OldData is stored (index ordered) index newList
        int [] indexNew =new int[newData.length];
        int [] indexOld =new int[sortedData.length];

        int index2=sortedData.length;
        double pivot2=newData[0];
        double pivot1=sortedData[0];
        int k1 = 0;

        for (int j = 0; j < newLength; j++){
            if (index2<newLength && k1<sortedData.length) {
                pivot1=sortedData[k1];
                pivot2=newData[index2-sortedData.length];

                if (pivot1<pivot2) {
                    newList[j]=pivot1;
                    newOrder[j]=k1;
                    indexOld[k1]=j;
                    k1++;
                } else {
                    newList[j]=pivot2;
                    newOrder[j]=index2;
                    indexNew[index2-sortedData.length]=j;
                    index2++;
                }

            } else if (index2<newLength) {
                newList[j]=newData[index2-sortedData.length];
                newOrder[j]=index2;
                indexNew[index2-sortedData.length]=j;
                index2++;
            } else {
                newList[j]=sortedData[k1];
                newOrder[j]=k1;
                indexOld[k1]=j;
                k1++;
            }
        }

        return new QuadupleGP(newList,newOrder,indexNew,indexOld);
    }

    //    Returns the index in the Neighbors data that contain new data
    private static PairIndex SubIndex (int [] Index, int numberOldData, int numberNewData) {
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

    private static class CSVstats {
        public ArrayList<double[]> changepoints, Gvalues;
        public double[] lambdas, precisions, tmrcas;

        public CSVstats(ArrayList<ArrayList<Double>> changepoints,
                        ArrayList<ArrayList<Double>> Gvalues,
                        ArrayList<Double> lambdas,
                        ArrayList<Double> precisions,
                        ArrayList<Double> tmrcas) {
            this.changepoints = new ArrayList<double[]>();
            for (ArrayList<Double> al : changepoints) {
                double[] arr = new double[al.size()];
                for (int i = 0; i < al.size(); i++) {
                    arr[i] = al.get(i);
                }
                this.changepoints.add(arr);
            }

            this.Gvalues = new ArrayList<double[]>();
            for (ArrayList<Double> al : Gvalues) {
                double[] arr = new double[al.size()];
                for (int i = 0; i < al.size(); i++) {
                    arr[i] = al.get(i);
                }
                this.Gvalues.add(arr);
            }

            this.lambdas = new double[lambdas.size()];
            for (int i = 0; i < lambdas.size(); i++) {
                this.lambdas[i] = lambdas.get(i);
            }

            this.precisions = new double[precisions.size()];
            for (int i = 0; i < precisions.size(); i++) {
                this.precisions[i] = precisions.get(i);
            }

            this.tmrcas = new double[tmrcas.size()];
            for (int i = 0; i < tmrcas.size(); i++) {
                this.tmrcas[i] = tmrcas.get(i);
            }
        }

        public int getSize() {
            int result = -1;
            if (changepoints.size() == Gvalues.size() && changepoints.size() == precisions.length &&
                    changepoints.size() == tmrcas.length) {
                result = changepoints.size();
            }
            return result;
        }
    }

    private static class TripGP {
        private double [] data;
        private int[] order;
        private int[] newOrder;

        public TripGP (double [] data, int [] order, int[] newOrder){
            this.data = data;
            this.order = order;
            this.newOrder = newOrder;
        }

        public double [] getData() {
            return data;
        }

        public int[] getOrder() {
            return order;
        }

        public int[] getNewOrder() {
            return newOrder;
        }
    }

    private static class QuadupleGP{
        private double[] data;
        private int[] order;
        private int[] positionNew;
        private int[] positionOld;

        public QuadupleGP (double [] data, int [] order, int [] positionNew, int[] positionOld){
            this.data = data;
            this.order = order;
            this.positionNew = positionNew;
            this.positionOld = positionOld;

        }
        public double[] getData() {return data;}
        public int[] getOrder() {return order;}
        public int [] getPositionNew() {return positionNew;}
        public int [] getPositionOld() {return positionOld;}
    }

    private static class PairIndex{
        private int[] orderNew;
        private int[] orderOld;

        public PairIndex (int [] orderNew, int [] orderOld){
            this.orderNew = orderNew;
            this.orderOld = orderOld;

        }
        public int[] getOrderNew() {return orderNew;}
        public int[] getOrderOld() {return orderOld;}
    }
}
