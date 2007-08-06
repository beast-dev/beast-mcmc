package dr.evomodel.continuous;

import no.uib.cipr.matrix.*;

import java.io.*;
import java.util.ArrayList;
import java.util.List;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;

/**
 * @author Marc Suchard
 */
public class BivariateDiscreteMapMaker {

    public BivariateDiscreteMapMaker(double[][] map) {
        xDim = map.length;
        yDim = map[0].length;
        this.map = map;
        totalDim = xDim * yDim;
        setUpRandomWalkMatrix();
    }

    private void setUpRandomWalkMatrix() {

        List<GraphEdge> edgeList = new ArrayList<GraphEdge>();
        for (int x = 0; x < xDim; x++) {
            for (int y = 0; y < yDim; y++) {
                int thisVertex = getIndex(x, y);
                //edgeList.add(new GraphEdge())
//					System.err.println("X = "+x+" y = "+y+" index = "+index);
                if (x < xDim - 1)
                    edgeList.add(new GraphEdge(
                            thisVertex, getIndex(x + 1, y), getWeight(x, y, x + 1, y)
                    ));
                if (x > 0)
                    edgeList.add(new GraphEdge(
                            thisVertex, getIndex(x - 1, y), getWeight(x, y, x - 1, y)
                    ));
                if (y < yDim - 1)
                    edgeList.add(new GraphEdge(
                            thisVertex, getIndex(x, y + 1), getWeight(x, y, x, y + 1)
                    ));
                if (y > 0)
                    edgeList.add(new GraphEdge(
                            thisVertex, getIndex(x, y - 1), getWeight(x, y, x, y - 1)
                    ));
            }
        }

        bandMatrix = new UpperSymmBandMatrix(totalDim, Math.max(xDim, yDim));
//				1);
        double[] rowTotal = new double[totalDim];
        randomWalkMatrix = new double[totalDim][totalDim];
        for (GraphEdge edge : edgeList) {
            int start = edge.getStart();
            int stop = edge.getStop();
            double weight = edge.getWeight();
            randomWalkMatrix[start][stop] = weight;
            rowTotal[start] -= weight;
            if (start < stop)
                bandMatrix.set(start, stop, weight);
        }
        for (int i = 0; i < totalDim; i++) {
            randomWalkMatrix[i][i] = rowTotal[i];
            bandMatrix.set(i, i, rowTotal[i]);
        }

//		System.err.println(new dr.math.matrixAlgebra.Matrix(randomWalkMatrix).toString());
//		System.err.println(bandMatrix.toString());
//		System.exit(-1);

    }


    public void doStuff() {

        SymmBandEVD evd = null;
        try {
            evd = SymmBandEVD.factorize(bandMatrix, totalDim);
        } catch (NotConvergedException e) {
            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
        }

        eVal = evd.getEigenvalues();
        DenseMatrix eigenvectors = evd.getEigenvectors();

//		System.err.println(bandMatrix.toString());
        eVec = new double[totalDim][totalDim];
        for (int i = 0; i < totalDim; i++) {
            for (int j = 0; j < totalDim; j++)
                eVec[i][j] = eigenvectors.get(i, j);
        }

        System.err.println("00 = " + getProbability(0, 0, 1000));
        System.err.println("01 = " + getProbability(0, 1, 1000));


        try {
            writeEigenvalues("eigenvalues.mat.gz", eVal);
            writeEigenvectors("eigenvectors.mat.gz", eVec);

            double[] newVec = readEigenvalues("eigenvalues.mat.gz");
            double[][] newMat = readEigenvectors("eigenvectors.mat.gz");

            System.err.println("E1 = " + newVec[0]);
            System.err.println("MAT11 = " + newMat[0][0]);
        } catch (IOException e) {
            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
        }

//		System.err.println("");

        System.exit(-1);


        System.err.println(eigenvectors.toString());

    }

    public static void writeEigenvectors(String file, double[][] mat) throws IOException {

        PrintWriter writer;
        if (file.endsWith("gz"))
            writer = new PrintWriter(
                    new OutputStreamWriter
                            (new GZIPOutputStream
                                    (new FileOutputStream(file))));
        else
            writer = new PrintWriter(new FileWriter(file));

        final int dim = mat.length;
        writer.println(dim);
        for (int i = 0; i < dim; i++) {
            for (int j = 0; j < dim; j++)
                writer.println(mat[i][j]);
        }
//		for(double d : vec )
//			writer.println(d);
        writer.close();

    }

    public static void writeEigenvalues(String file, double[] vec) throws IOException {

        PrintWriter writer;
        if (file.endsWith("gz"))
            writer = new PrintWriter(
                    new OutputStreamWriter
                            (new GZIPOutputStream
                                    (new FileOutputStream(file))));
        else
            writer = new PrintWriter(new FileWriter(file));

        writer.println(vec.length);
        for (double d : vec)
            writer.println(d);
        writer.close();

    }

    public static double[] readEigenvalues(String file) throws IOException {

        BufferedReader reader;
        if (file.endsWith("gz"))
            reader = new BufferedReader
                    (new InputStreamReader
                            (new GZIPInputStream
                                    (new FileInputStream(file))));
        else
            reader = new BufferedReader(new FileReader(file));

        final int dim = Integer.parseInt(reader.readLine());
        double[] vec = new double[dim];
        for (int i = 0; i < dim; i++)
            vec[i] = Double.parseDouble(reader.readLine());
        reader.close();

        return vec;

    }

    public static double[][] readEigenvectors(String file) throws IOException {

        BufferedReader reader;
        if (file.endsWith("gz"))
            reader = new BufferedReader
                    (new InputStreamReader
                            (new GZIPInputStream
                                    (new FileInputStream(file))));
        else
            reader = new BufferedReader(new FileReader(file));

        final int dim = Integer.parseInt(reader.readLine());
        double[][] mat = new double[dim][dim];
        for (int i = 0; i < dim; i++) {
            for (int j = 0; j < dim; j++)
                mat[i][j] = Double.parseDouble(reader.readLine());
        }
        reader.close();

        return mat;

    }


    private double getProbability(int I, int J, double time) {
        double probability = 0;
        for (int k = 0; k < totalDim; k++) {
            probability += eVec[I][k] * Math.exp(time * eVal[k]) * eVec[J][k];
        }
        return probability;

    }


    public class GraphEdge {
        private int start;
        private int stop;
        private double weight;

        public GraphEdge(int start, int stop, double weight) {
            this.start = start;
            this.stop = stop;
            this.weight = weight;
        }

        public int getStart() {
            return start;
        }

        public void setStart(int start) {
            this.start = start;
        }

        public int getStop() {
            return stop;
        }

        public void setStop(int stop) {
            this.stop = stop;
        }

        public double getWeight() {
            return weight;
        }

        public void setWeight(double weight) {
            this.weight = weight;
        }
    }

    public int getIndex(int x, int y) {
        return x * yDim + y;
    }

    public double getWeight(int x0, int y0, int x1, int y1) {
        return 1.0 / (100 * Math.abs(map[x0][y0] - map[x1][y1]) + 1.0);
    }


    public static void main(String[] args) {
        new BivariateDiscreteMapMaker(sillyMap).doStuff();


    }

    private double[][] map;
    private double[][] randomWalkMatrix;
    private double[] eVal;
    private double[][] eVec;
    Matrix bandMatrix;
    private int xDim;
    private int yDim;
    private int totalDim;

    public static double[][] sillyMap = {
            {0.0, 100.0, 0.0, 100.0, 0.0},
            {0.0, 0.0, 0.0, 100.0, 0.0},
            {100.0, 0.0, 100.0, 100.0, 0.0},
            {100.0, 0.0, 100.0, 100.0, 0.0},
            {100.0, 0.0, 0.0, 0.0, 0.0}
    };

}
