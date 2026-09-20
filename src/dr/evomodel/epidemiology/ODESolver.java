package dr.evomodel.epidemiology;

import org.apache.commons.math.ode.DerivativeException;
import org.apache.commons.math.ode.FirstOrderDifferentialEquations;
import org.apache.commons.math.ode.FirstOrderIntegrator;
import org.apache.commons.math.ode.IntegratorException;
import org.apache.commons.math.ode.events.EventHandler;
import org.apache.commons.math.ode.nonstiff.DormandPrince853Integrator;
import org.apache.commons.math.ode.sampling.FixedStepHandler;
import org.apache.commons.math.ode.sampling.StepNormalizer;

public class ODESolver extends CompartmentalModelSimulator {

    private final double absTolerance;
    private final double relTolerance;

    public ODESolver(CompartmentalModel compartmentalModel) {
        super(compartmentalModel);
        this.absTolerance = 1e-6;
        this.relTolerance = 1e-6;
    }

    public ODESolver(CompartmentalModel compartmentalModel,
                     double absTolerance,
                     double relTolerance) {
        super(compartmentalModel);
        this.absTolerance = absTolerance;
        this.relTolerance = relTolerance;
    }

    // Solves ODE system and records compartment counts at each grid point
    // ODE integration proceeds in forward time, starting at
    // cutOff - oldestOrigin (both of which are in backward time), and goes
    // on until we reach present time, with forward time value numerically equally to cutOff
    @Override
    public void simulateTrajectory(){

        final SimulationState state = initializeSimulation();
        final double youngerForwardOrigTime = compartmentalModel.getYoungerForwardOrigTime();

        FirstOrderDifferentialEquations odes = new FirstOrderDifferentialEquations() {
            @Override
            public int getDimension() {
                return numSpecies;
            }

            @Override
            public void computeDerivatives(double t, double[] y, double[] yDot) throws DerivativeException {
                double[] derivatives = compartmentalModel.getCompartmentDerivatives(y, t);
                System.arraycopy(derivatives, 0, yDot, 0, numSpecies);
            }
        };

        // Integrator
        // Dormand-Prince RK8(5,3): adaptive step size, 13 function evaluations
        // per step, mixed absolute/relative error control
        FirstOrderIntegrator integrator = new DormandPrince853Integrator(
                1e-10, // minimum step size
                intervalWidth, // maximum step size
                absTolerance,
                relTolerance);


        // handle second pathogen introduction
        // the introduction creates a discontinuity (one individual moves from SS to SI or IS)
        // this locates the time exactly and restarts the integrator from the modified state
        if(!compartmentalModel.isSecondPathogenIntroduced()) {
            integrator.addEventHandler(new EventHandler() {
                //@Override
                //public void init(double t0, double[] y0, double t) {}

                @Override
                public double g(double t, double[] y) {
                    // zero-crossing at s == youngerForwardOrigTime
                    return t - youngerForwardOrigTime;
                }

                @Override
                public int eventOccurred(double t, double[] y, boolean increasing) {
                    // introduce second pathogen by modifying y
                    // SS reduced by 1, IS or SI increased by 1
                    compartmentalModel.introduceSecondPathogen(t, y);

                    // keep state.currentCounts consistent with y so that
                    // recordCompartmentCountsUpTo uses the correct counts
                    System.arraycopy(y, 0, state.currentCounts, 0, numSpecies);

                    // RESET_STATE causes the integrator to restart from the modified y
                    return EventHandler.RESET_STATE;
                }

                @Override
                public void resetState(double t, double[] y) {
                    // y has already been modified in eventOccurred
                    // do nothing
                }
                // max interval between event checks, root-find convergence tolerance,
                // and max root-finding iterations
            }, intervalWidth, 1e-10, 100);
        }

        // step handler to record compartment counts at grid points
        // StepNormalizer calls handleStep at regular intervals of width
        // intervalWidth/100 so that we get less discretization than
        // grid spacing would impose
        integrator.addStepHandler(new StepNormalizer(
                intervalWidth / 100.0,
                new FixedStepHandler() {
                    @Override
                    public void handleStep(double t, double[] y, double[] yDot,
                                           boolean isLast) throws DerivativeException {
                        System.arraycopy(y, 0, state.currentCounts, 0, numSpecies);
                        state.simulationTime = t;
                        recordCompartmentCountsUpTo(state, t);
                    }
                }
        ));

        // run integrator from oldest origin to present time
        double[] y = state.currentCounts.clone();
        // cutOff - oldestOrigin
        double t0 = state.simulationTime;
        // end of simulation
        double t1 = cutOff;

        try{
            integrator.integrate(odes, t0, y, t1, y);
        } catch (IntegratorException e){
            // integration failed
            // set lineageConstraintViolated = true;
            // will cause proposal to be rejected
            lineageConstraintViolated = true;
        } catch (DerivativeException e){
            // derivative computation failed
            lineageConstraintViolated = true;
        }
    }
}
