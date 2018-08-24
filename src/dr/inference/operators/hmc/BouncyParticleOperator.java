package dr.inference.operators.hmc;

import dr.inference.hmc.GradientWrtParameterProvider;
import dr.inference.model.Likelihood;
import dr.inference.model.Parameter;
import dr.inference.operators.*;
import dr.math.MathUtils;
import dr.math.distributions.NormalDistribution;
import dr.util.Transform;

import java.util.Arrays;
import java.io.FileWriter;
import java.io.FileOutputStream;
import java.io.IOException;

import static org.netlib.lapack.Dlacon.j;

public class BouncyParticleOperator {

    private double t; //todo randomize the length a little bit.
    final double t0;
    private double[] phi_w;
    private double[] mu;
    private double[] location;
    private double[] v;
    private double[][] precisionmatrix;
    final NormalDistribution drawDistribution;

    public static void main(String[] args)throws Exception{

        FileWriter fw = new FileWriter("data.csv");
        int iterations = 1000000;
        int dim = 10;
        double[][] result = new double[iterations][dim];
        BouncyParticleOperator BPS1 = new BouncyParticleOperator();

        for (int i = 0; i < iterations; i ++){

            BPS1.bpsOneStep();
            BPS1.t = BPS1.t0;
            for (int j = 0; j < dim; j ++){
                result[i][j] = BPS1.location[j];
            }

        }



        try{
        for (int i = 0; i < iterations; i ++) {
            for (int j = 0; j < dim; j ++) {
                fw.write(result[i][j] + ",");
            }
            fw.write("\n");
        }
        fw.close();

        } catch (IOException e) {
            e.printStackTrace();
        }

    }


 //todo resume the constructor
//    public BouncyParticleOperator (double[][] precisionmatrix,double[] location, double[] mu, double t) {
//
//        this.location = location;
//        this.mu = mu;
//        this.t = t;
//        this.precisionmatrix = precisionmatrix;
//        this.drawDistribution = new NormalDistribution(0, Math.sqrt(1));
//        this.v = drawV(drawDistribution, location.length);
//    }

        public BouncyParticleOperator () {


            double[][] precisionmatrix = new double[][] {
                    {11.1016047,-12.0787204, 5.58064211,-3.1689313, 1.094961,-4.961749,5.11981149,-1.5834688, 3.5623476,-0.7376637},
                    { -12.0787204,18.3400620,-11.73593592, 0.9061316,-2.016980, 6.451601,-1.05863904, 0.2369949,-3.8070834,-0.5602705},
                    { 5.5806421,-11.7359359,28.49214394,-20.7446633, 8.302015,-4.912537,-0.03486475, 3.2872308, 3.4449720,-2.2412823},
                    {-3.1689313, 0.9061316,-20.74466331,40.7345307,-18.290551, 2.218922,-7.73442938,-4.2918288,-0.3767421,8.7189560},
                    { 1.0949613,-2.0169801, 8.30201463,-18.2905513,21.939836,-10.470282,-3.05356573, 1.8893231, 5.0052153,-2.9950543},
                    {-4.9617490, 6.4516010,-4.91253689, 2.2189222,-10.470282,20.157010,-7.83388716, 5.2679596,-12.3805480,1.4309123},
                    { 5.1198115,-1.0586390,-0.03486475,-7.7344294,-3.053566,-7.833887, 25.75186700,-9.0823931, 4.4477436,-4.5515479},
                    {-1.5834688, 0.2369949, 3.28723080,-4.2918288, 1.889323, 5.267960,-9.08239305,13.3924830,-12.0389557,0.3323174},
                    { 3.5623476,-3.8070834, 3.44497199,-0.3767421, 5.005215,-12.380548,4.44774359,-12.0389557,25.0322293,-6.4872396},
                    {-0.7376637,-0.5602705,-2.24128232, 8.7189560,-2.995054, 1.430912,-4.55154794, 0.3323174,-6.4872396,8.1990433}
            };
            double[] location = new double[]{1.0973772,1.0654628,0.9832338,0.9924755,1.1105650,1.1999598,0.8974071, 1.0889359,0.8875886,1.0883105};
            double[] mu = new double[]{1.0460811,1.1343102,1.0997168,0.9632376,0.9868017,1.1777099,0.9535731,0.9826153,0.9595986,1.1414815};

            double t = 2.131453;
            double t0 = 2.131453;

            this.location = location;
            this.mu = mu;
            this.t = t;
            this.t0 = t0;
            this.precisionmatrix = precisionmatrix;
            this.drawDistribution = new NormalDistribution(0, Math.sqrt(1));

    }

    public double doOperation(double[] location, double[] v, double t, double[] phi_w) {


        bpsOneStep();

        return 0.0;
    }

    public void bpsOneStep() {

        v = drawV(drawDistribution, location.length);

        phi_w = matrixMultiplier(precisionmatrix, addArray(location, mu, true));
        while (t > 0) {

            double[] w = addArray(location, mu, true);
            double[] phi_v = matrixMultiplier(precisionmatrix, v);
            double w_phi_w = getDotProduct(w, phi_w);//todo multiple precision matrix should be a class.
            double v_phi_w = getDotProduct(v, phi_w);//todo use a construct to store all of the temporary values.
            double v_phi_v = getDotProduct(v, phi_v);

            double tMin = Math.max(0.0, - v_phi_w/v_phi_v);

            double U_min = energyProvider(phi_w, phi_v, w, tMin);


            if( Double.isNaN(v_phi_w)){

                System.exit(-99);
            }
            double bounceTime = getBounceTime(v_phi_v, v_phi_w, U_min, w_phi_w);


                    TravelTime time_to_bdry = getTimeToBoundary(location, v);


            bpsUpdate(time_to_bdry, bounceTime, t);



        }


    }

    public void bpsUpdate(TravelTime traveltime, double bouncetime, double timeRemain) {

        if (timeRemain < Math.min(traveltime.minTime, bouncetime)) { //no event

            location = addArray(location, getConstantVector(v, timeRemain), false);
            t = 0.0;


        } else if (traveltime.minTime <= bouncetime) { //against the boundary

            location = addArray(location, getConstantVector(v, traveltime.minTime), false);
            location[traveltime.minIndex] = 0.0;
            v[traveltime.minIndex] = getConstantVector(v, -1.0)[traveltime.minIndex];
            t -=  traveltime.minTime;
            phi_w = matrixMultiplier(precisionmatrix, addArray(location, mu, true)); //todo to deal with the duplicated code for updating phi_w


        } else { //against the gradient

            location = addArray(location, getConstantVector(v, bouncetime), false);
            phi_w = matrixMultiplier(precisionmatrix,addArray(location, mu, true));
            v = bounceAgainst(v, phi_w);
            t -= bouncetime;

        }

    }

    public double[] bounceAgainst(double[] v, double[] minusGrad) {

            double[] finalV;
            double[] verticalV;

            verticalV = getConstantVector(minusGrad, getDotProduct(v, minusGrad)/getDotProduct(minusGrad, minusGrad));
            finalV = addArray(v, getConstantVector(verticalV, -2.0), false);

            return finalV;

    }

    public double getBounceTime(double v_phi_v, double v_phi_w, double U_min, double w_phi_w) {

        double a = v_phi_v;
        double b = 2 * v_phi_w;
        double c = 2.0 * Math.log(1 - MathUtils.nextDouble()) - U_min + w_phi_w; //TODO CHCEK IF ROOT EXISTS

        return (- b + Math.sqrt(b * b - 4 * a * c))/2/a;

    }

    public double energyProvider(double[] phi_w, double[] phi_v, double[] w, double t){

        return getDotProduct(addArray(w, getConstantVector(v, t), false),
                addArray(phi_w, getConstantVector(phi_v, t), false));

    }

    static double[] drawV(final NormalDistribution distribution, final int dim) {

        double[] v = new double[dim];
        for (int i = 0; i < dim; i++) {
            v[i] = (Double) distribution.nextRandom();
        }
        return v;
    }

    public double[] matrixMultiplier(double[][] A, double[] B){

        final int dim = B.length;
        double[] mResult = new double[dim];

        for (int i = 0; i < dim; i ++){
            for (int j = 0; j< dim; j ++) {
                mResult[i] += A[i][j] * B[j];
            }

        }
        return mResult;
    }

    public TravelTime getTimeToBoundary(double[] location, double[] v) {

        double[] travelTime = elementWiseVectorDivisionAbsValue(location, v);

        int index = 0;

        double minTime = Double.MAX_VALUE;

        for (int i = 0; i < travelTime.length; i ++) {
            if (travelTime[i] != 0.0 && location[i] * v[i] < 0) {
                if (travelTime[i] < minTime) {
                   index = i;
                   minTime = travelTime[i];
                }
            }
        }

        return new TravelTime(travelTime, minTime, index);

    }

    public double[] elementWiseVectorDivisionAbsValue(double[] location, double[] y) {

        final int dim = location.length;
        double[] z = new double[dim];
        for (int i = 0; i < dim; i++) {
            z[i] = Math.abs(location[i] / y[i]); //todo make sure b[i] != 0
        }
        return z;
    }

    public static double[] addArray(double[] a, double[] b, boolean subtract) {

        assert (a.length == b.length);
        final int dim = a.length;

        double result[] = new double[dim];
        for (int i = 0; i < dim; i++) {
            if(!subtract) {
                result[i] = a[i] + b[i];
            } else {
                result[i] = a[i] - b[i];
            }

        }

        return result;
    }

    public static double getDotProduct(double[] a, double[] b) { //todo get rid of code duplication

        assert (a.length == b.length);
        final int dim = a.length;

        double total = 0.0;
        for (int i = 0; i < dim; i++) {

            total += a[i]*b[i];
        }
        return total;
    }

    public static double[] getConstantVector(double[] a, double c){

        double[] cx = new double[a.length];
        for (int i = 0; i < a.length; i ++){
            cx[i] = a[i] * c;
        }

        return cx;
    }

    public class TravelTime {

        double[] traveltime;
        double minTime;
        int minIndex;

        public TravelTime (double[] traveltime, double minTime, int minIndex){
            this.traveltime = traveltime;
            this.minTime = minTime;
            this.minIndex = minIndex;
        }
    }

}
