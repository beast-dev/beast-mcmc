/*
 * NewHamiltonianMonteCarloOperator.java
 *
 * Copyright (c) 2002-2017 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.inference.operators.hmc;

import dr.inference.hmc.GradientWrtParameterProvider;
import dr.inference.hmc.PrecisionMatrixVectorProductProvider;
import dr.inference.model.Parameter;
import dr.inference.operators.GibbsOperator;
import dr.inference.operators.SimpleMCMCOperator;
import dr.math.MathUtils;
import dr.math.distributions.NormalDistribution;
import dr.math.matrixAlgebra.ReadableVector;
import dr.math.matrixAlgebra.WrappedVector;

import java.util.Arrays;

/**
 * @author Zhenyu Zhang
 * @author Marc A. Suchard
 */

public class NewBouncyParticleOperator extends SimpleMCMCOperator implements GibbsOperator {

    private final NormalDistribution drawDistribution;
    private final Parameter parameter;

    @Override
    public String getPerformanceSuggestion() {
        return null;
    }

    @Override
    public String getOperatorName() {
        return "Bouncy Particle operator";
    }

    private static final boolean DEBUG = false;

    @Override
    public double doOperation() {
        bpsOneStep();
        return 0.0;
    }

    private void bpsOneStep() {

        WrappedVector position = new WrappedVector.Raw(getInitialPosition());

        WrappedVector velocity = new WrappedVector.Raw(drawVelocity(drawDistribution, masses));

        WrappedVector negativeGradient = getNegativeGradient(position);

        double remainingTime = travelTime;
        while (remainingTime > 0) {

            ReadableVector Phi_v = getPrecisionProduct(velocity);

//            double[] mu = getMU();
//            double[] w = addArray(location, mu, true);
//            double[] phi_v = matrixMultiplier(precisionMatrix, velocity);

            double v_Phi_x = innerProduct(velocity, negativeGradient);
            double v_Phi_v = innerProduct(velocity, Phi_v);



//            double w_phi_w = getDotProduct(w, phi_w);//todo multiple precision matrix should be a class.
//            double v_phi_w = getDotProduct(velocity, phi_w);//todo use a construct to store all of the temporary values.
//            double v_phi_v = getDotProduct(velocity, phi_v);

            double tMin = Math.max(0.0, - v_Phi_x / v_Phi_v);
            double U_min = tMin * tMin / 2 * v_Phi_v + tMin * v_Phi_x;
//            double U_min = energyProvider(velocity, phi_w, phi_v, w, tMin);

//            if( Double.isNaN(v_phi_w)){
//
//                System.exit(-99);
//            }

            double bounceTime = getBounceTime(v_Phi_v, v_Phi_x, U_min);
            TravelTime time_to_bdry = getTimeToBoundary(position, velocity);

            remainingTime = doBounce(
                    remainingTime, bounceTime, time_to_bdry,
                    position, velocity, negativeGradient, Phi_v
            );

            //System.err.println("location is (inside) +  " + Arrays.toString(location));

        }

        setParameter(position);

        //System.err.println("location is (outside) +  " + Arrays.toString(location));
    }

    private double doBounce(double remainingTime, double bounceTime, TravelTime timeToBoundary,
                            WrappedVector position, WrappedVector velocity,
                            WrappedVector negativeGradient, ReadableVector Phi_v) {

        if (remainingTime < Math.min(timeToBoundary.minTime, bounceTime)) { // No event during remaining time

            updatePosition(position, velocity, remainingTime);
            remainingTime = 0.0;


        } else if (timeToBoundary.minTime < bounceTime) { // Bounce against the boundary

            updatePosition(position, velocity, timeToBoundary.minTime);
            updateNegativeGradient(negativeGradient, timeToBoundary.minTime, Phi_v);

            position.set(timeToBoundary.minIndex, 0.0);
            velocity.set(timeToBoundary.minIndex, -1 * velocity.get(timeToBoundary.minIndex));

            remainingTime -= timeToBoundary.minTime;

        } else { // Bounce caused by the gradient

            updatePosition(position, velocity, bounceTime);
            updateNegativeGradient(negativeGradient, bounceTime, Phi_v);
            updateVelocity(velocity, negativeGradient, new WrappedVector.Raw(masses));

            remainingTime -= bounceTime;

        }

        return remainingTime;
    }

    private void setParameter(ReadableVector position) {
        for (int j = 0, dim = position.getDim(); j < dim; ++j) {
            parameter.setParameterValueQuietly(j, position.get(j));
        }
        parameter.fireParameterChangedEvent();
    }

    private void updateVelocity(WrappedVector velocity, WrappedVector negativeGradient, ReadableVector masses) {
        // TODO Handle masses

        double vg = innerProduct(velocity, negativeGradient); // TODO Isn't this already computed
        double gg = innerProduct(negativeGradient, negativeGradient);

        for (int i = 0, len = velocity.getDim(); i < len; ++i) {
            velocity.set(i,
                    velocity.get(i) + 2 * vg / gg * negativeGradient.get(i));
        }
    }

    private void updateNegativeGradient(WrappedVector negativeGradient, double time, ReadableVector Phi_v) {
        for (int i = 0, len = negativeGradient.getDim(); i < len; ++i) {
            negativeGradient.set(i, negativeGradient.get(i) + time * Phi_v.get(i));
        }
    }

    private void updatePosition(WrappedVector position, WrappedVector velocity, double time) {
        for (int i = 0, len = position.getDim(); i < len; ++i) {
            position.set(i, position.get(i) + time * velocity.get(i));
        }
    }
    
    private double getBounceTime(double v_phi_v, double v_phi_x, double u_min) {
        double a = v_phi_v / 2;
        double b = v_phi_x;
        double c = MathUtils.nextExponential(1) - u_min;
        return (-b + Math.sqrt(b * b - 4 * a * c));
    }

    private WrappedVector getNegativeGradient(ReadableVector position) {

        setParameter(position);

        double[] gradient = gradientProvider.getGradientLogDensity();
        for (int i = 0, len = gradient.length; i < len; ++i) {
            gradient[i] = -1 * gradient[i];
        }

        return new WrappedVector.Raw(gradient);
    }

    private double innerProduct(ReadableVector x, ReadableVector y) {

        assert (x.getDim() == y.getDim());

        double sum = 0;
        for (int i = 0, dim = x.getDim(); i < dim; ++i) {
            sum += x.get(i) * y.get(i);
        }

        return sum;
    }

    private ReadableVector getPrecisionProduct(ReadableVector velocity) {

        setParameter(velocity);

        double[] product = multiplicationProvider.getProduct(parameter);

        return new WrappedVector.Raw(product);
    }

    private static double[] drawVelocity(final NormalDistribution distribution, double[] masses) {

        double[] velocity = new double[masses.length];

        for (int i = 0; i < masses.length; i++) {
            velocity[i] = (Double) distribution.nextRandom() / Math.sqrt(masses[i]);
        }
        return velocity;
    }

    private TravelTime getTimeToBoundary(ReadableVector position, ReadableVector velocity) {

        assert (position.getDim() == velocity.getDim());

        int index = -1;
        double minTime = Double.MAX_VALUE;

        for (int i = 0, len = position.getDim(); i < len; ++i) {

            double travelTime = Math.abs(position.get(i) / velocity.get(i));

            if (travelTime > 0.0 && headingAwayFromBoundary(position.get(i), velocity.get(i))) {

                if (travelTime < minTime) {
                    index = i;
                    minTime = travelTime;
                }
            }
        }

        return new TravelTime(minTime, index);

    }

    private boolean headingAwayFromBoundary(double position, double velocity) {
        return position * velocity < 0.0;
    }

    private class TravelTime {

//        double[] traveltime;
        double minTime;
        int minIndex;

        private TravelTime (//double[] traveltime,
                           double minTime, int minIndex){
//            this.traveltime = traveltime;
            this.minTime = minTime;
            this.minIndex = minIndex;
        }
    }

    private double[] getInitialPosition() {
        return parameter.getParameterValues();
    }

    public NewBouncyParticleOperator(GradientWrtParameterProvider gradientProvider,
                                     PrecisionMatrixVectorProductProvider multiplicationProvider,
                                     double weight) {

        this.gradientProvider = gradientProvider;
        this.multiplicationProvider = multiplicationProvider;
        this.parameter = gradientProvider.getParameter();
        this.drawDistribution = new NormalDistribution(0, 1);

        setWeight(weight);
        checkParameterBounds(parameter);

        masses = setupPreconditionedMatrix();

        // TODO Determine travelTime
        travelTime = 0.05;
    }

    private static void checkParameterBounds(Parameter parameter) {
        // TODO
    }

    private double[] setupPreconditionedMatrix() {
        double[] masses = new double[parameter.getDimension()];
        Arrays.fill(masses, 1.0); // TODO
        return masses;
    }

    private double travelTime;

    private final double[] masses;

    private GradientWrtParameterProvider gradientProvider;
    private PrecisionMatrixVectorProductProvider multiplicationProvider;
}
