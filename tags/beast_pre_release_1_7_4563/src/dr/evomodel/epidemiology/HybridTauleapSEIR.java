package dr.evomodel.epidemiology;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;
import cern.jet.random.Poisson;

/**
 * Class implementing Cao et al.'s hybrid tau-leaping algorithm
 *
 * @author tvaughan
 * @author Trevor Bedford
 */
public class HybridTauleapSEIR {

	SEIRState initialState;
	double exposeRate, infectRate, recoverRate;
	double alpha;
	Random random;

	boolean useExposed;

	/**
	 * Constructor
	 *
	 * @param state0	initial state of system
	 * @param infect	infection rate
	 * @param recover	recovery rate
	 * @param alpha		threshold for determining critical reactions
	 */
	public HybridTauleapSEIR(SEIRState state0,
			double expose, double infect, double recover,
			boolean useExposed,
			double alpha, Random random) {
		super();
		this.initialState = state0.copy();
		this.exposeRate = expose;
		this.infectRate = infect;
		this.recoverRate = recover;
		this.useExposed = useExposed;
		this.alpha = alpha;
		this.random = random;
	}

	/**
	 * Set system state
	 *
	 * @param newState
	 */
	public void setState(SEIRState newState) {
		initialState = newState.copy();
	}

	/**
	 * Perform one time-step of fixed length
	 * (Does not use exposed compartment.)
	 *
	 * @param	dt	time-step size
	 */
	public boolean step(double dt) {

		double t = 0.0;

		boolean critical_step = false;

		while (true) {

			// Calculate propensities:
			double a_infect = infectRate*initialState.S*initialState.I;
			double a_recover = recoverRate*initialState.I;

			// Determine which reactions are "critical"
			// and calculate next reaction time:
			double dtCrit = dt - t;
			int critReactionToFire = 0;

			boolean infectIsCrit = false;
			if (initialState.S < a_infect*dtCrit + alpha*Math.sqrt(a_infect*dtCrit)) {

				// Infection reaction is critical
				infectIsCrit = true;

				// Determine whether infection will fire:
				double thisdt = -Math.log(random.nextDouble())/a_infect;
				if (thisdt < dtCrit) {
					dtCrit = thisdt;
					critReactionToFire = 1;
				}
			}

			boolean recoverIsCrit = false;
			if (initialState.I < a_recover*dtCrit + alpha*Math.sqrt(a_recover*dtCrit)) {

				// Recovery is critical:
				recoverIsCrit = true;

				// Determine whether recovery will fire:
				double thisdt = -Math.log(random.nextDouble())/a_recover;
				if (thisdt < dtCrit) {
					dtCrit = thisdt;
					critReactionToFire = 2;
				}
			}

			// Reaction has been marked as critical, so this is a
                        // critical step:
                        if (infectIsCrit || recoverIsCrit)
                            critical_step = true;

			// Update time:
			t += dtCrit;

			// tau-leap non-critical reactions:
			if (infectIsCrit == false) {
				int q = Poisson.staticNextInt(a_infect*dtCrit);
				initialState.S -= q;
				initialState.I += q;
			}

			if (recoverIsCrit == false) {
				int q = Poisson.staticNextInt(a_recover*dtCrit);
				initialState.I -= q;
				initialState.R += q;
			}

			// End step if no critical reaction fires:
			if (critReactionToFire == 0)
				break;

			// implement one critical reaction:
			if (critReactionToFire == 1) {
				// Infection
				initialState.S -= 1;
				initialState.I += 1;
			} else {
				// Recovery
				initialState.I -= 1;
				initialState.R += 1;
			}

		}

		// Check for negative populations:
		if (initialState.S < 0 || initialState.I < 0 || initialState.R < 0) {
//			System.err.println("Error: negative population detected. Rejecting trajectory.");
            throw new RuntimeException("Error: negative population detected. Rejecting trajectory.");
//			System.exit(1);
		}


		// Update state time:
		initialState.time += dt;

		return critical_step;

	}

	/**
	 * Perform one fixed-size time step using exposed compartment
	 *
	 * @param dt	length of time step
	 * @return true or false depending on whether step involved "critical" reactions
	 */
	public boolean step_exposed(double dt) {
		double t = 0.0;

		boolean critical_step = false;

		while (true) {

			// Calculate propensities:
			double a_expose = exposeRate*initialState.S*initialState.I;
			double a_infect = infectRate*initialState.E;
			double a_recover = recoverRate*initialState.I;

			// Determine which reactions are "critical"
			// and calculate next reaction time:
			double dtCrit = dt - t;
			int critReactionToFire = 0;

			boolean exposeIsCrit = false;
			if (initialState.S < a_expose*dtCrit + alpha*Math.sqrt(a_expose*dtCrit)) {

				// Exposure reaction is critical
				exposeIsCrit = true;

				// Determine whether exposure will fire:
				double thisdt = -Math.log(random.nextDouble())/a_expose;
				if (thisdt < dtCrit) {
					dtCrit = thisdt;
					critReactionToFire = 1;
				}
			}

			boolean infectIsCrit = false;
			if (initialState.E < a_infect*dtCrit + alpha*Math.sqrt(a_infect*dtCrit)) {

				// Infection reaction is critical
				infectIsCrit = true;

				// Determine whether infection will fire:
				double thisdt = -Math.log(random.nextDouble())/a_infect;
				if (thisdt < dtCrit) {
					dtCrit = thisdt;
					critReactionToFire = 2;
				}
			}

			boolean recoverIsCrit = false;
			if (initialState.I < a_recover*dtCrit + alpha*Math.sqrt(a_recover*dtCrit)) {

				// Recovery is critical:
				recoverIsCrit = true;

				// Determine whether recovery will fire:
				double thisdt = -Math.log(random.nextDouble())/a_recover;
				if (thisdt < dtCrit) {
					dtCrit = thisdt;
					critReactionToFire = 3;
				}
			}

			// Reaction has been marked as critical, so this is a
                        // critical step:
                        if (exposeIsCrit || infectIsCrit || recoverIsCrit)
                            critical_step = true;

			// Update time:
			t += dtCrit;

			// tau-leap non-critical reactions:
			if (exposeIsCrit == false) {
				int q = Poisson.staticNextInt(a_expose*dtCrit);
				initialState.S -= q;
				initialState.E += q;
			}

			// tau-leap non-critical reactions:
			if (infectIsCrit == false) {
				int q = Poisson.staticNextInt(a_infect*dtCrit);
				initialState.E -= q;
				initialState.I += q;
			}

			if (recoverIsCrit == false) {
				int q = Poisson.staticNextInt(a_recover*dtCrit);
				initialState.I -= q;
				initialState.R += q;
			}

			// End step if no critical reaction fires:
			if (critReactionToFire == 0)
				break;

			// implement one critical reaction:
			switch (critReactionToFire) {
			case 1:
				// Exposure
				initialState.S -= 1;
				initialState.E += 1;
				break;

			case 2:
				// Infection
				initialState.E -= 1;
				initialState.I += 1;
				break;

			case 3:
				// Recovery
				initialState.I -= 1;
				initialState.R += 1;
				break;
			}

		}

		// Check for negative populations:
		if (initialState.S < 0 || initialState.E < 0 || initialState.I < 0 || initialState.R < 0) {
            System.err.println("Error: negative population detected. Rejecting trajectory.");
            throw new RuntimeException("Error: negative population detected. Rejecting trajectory.");
//			System.exit(1);
		}


		// Update state time:
		initialState.time += dt;

		return critical_step;

	}

	/**
	 * Generate trajectory
	 *
	 * @param T			Integration time
	 * @param Nt		Number of time-steps
	 * @param Nsamples	Number of samples to record
	 *
	 * @return List of SEIRState instances representing sampled trajectory.
	 */
	public List<SEIRState> genTrajectory(double T, int Nt, int Nsamples, List<Integer> criticalTrajectories) {

		// Determine time-step size:
		double dt = T/(Nt-1);

		// Determine number of time steps per sample:
		int stepsPerSample = (Nt-1)/(Nsamples-1);

		// Allocate memory for sampled states:
		List<SEIRState> trajectory = new ArrayList<SEIRState>();

		// Sample first state:
		trajectory.add(initialState.copy());

		for (int tidx=1; tidx<Nt; tidx++) {

			// Increment state:
			boolean crit = false;

			if (useExposed)
				crit = step_exposed(dt);
			else
				crit = step(dt);

			// Sample if necessary:
			if (tidx % stepsPerSample == 0) {
				trajectory.add(initialState.copy());

				// Adjust critical trajectory count:
//				if (crit)
//					criticalTrajectories.set(tidx/stepsPerSample,
//							criticalTrajectories.get(tidx/stepsPerSample)+1);
			}
		}

		return trajectory;

	}

	/**
	 * Generate trajectory and store trajectory in class
     * Store one state per time step
	 *
	 * @param T			Integration time
	 * @param Nt		Number of time-steps
	 *
	 */

	public void storeTrajectory(double T, int Nt) {

		// Determine time-step size:
		double dt = T/(Nt-1);

		// Allocate memory for sampled states:
		List<SEIRState> trajectory = new ArrayList<SEIRState>();

		// Sample first state:
		trajectory.add(initialState.copy());

		for (int tidx=1; tidx<Nt; tidx++) {

			// Increment state:
			boolean crit = false;

			if (useExposed)
				crit = step_exposed(dt);
			else
				crit = step(dt);

			// Sample
			trajectory.add(initialState.copy());
		}

	}


    // slightly modified method for BEAST Operator
    	public List<SEIRState> genTrajectory(double T, int Nt, int Nsamples, int ntaxa, Boolean check) {

		// Determine time-step size:
		double dt = T/(Nt-1);

		// Determine number of time steps per sample:
		int stepsPerSample = (Nt-1)/(Nsamples-1);

		// Allocate memory for sampled states:
		List<SEIRState> trajectory = new ArrayList<SEIRState>();

		// Sample first state:
		trajectory.add(initialState.copy());

		for (int tidx=1; tidx<Nt; tidx++) {


			if (useExposed)
				step_exposed(dt);
			else
				step(dt);

			// Sample if necessary:
			if (tidx % stepsPerSample == 0) {
                if ((initialState.I < 1) && (!useExposed || (initialState.E == 0)))
//                if (trajectory.get(trajectory.size()-1).S == state.S ) zeroCount++;
//                if (check && zeroCount > Nsamples/2.)
                    throw new RuntimeException("Abort simulation. No infecteds left.");

				trajectory.add(initialState.copy());

			}
		}

        if ((initialState.I < 1) && (!useExposed || (initialState.E == 0)))
            throw new RuntimeException("Abort simulation. No infecteds left.");

		return trajectory;

	}

	/**
	 * Main method: for debugging only
	 *
	 * @param args
	 */
	public static void main(String[] args) {


          double gamma = 1.08e-3; //.25;       // infected recovery rate
          double beta = 8.14e-3;//4.48E-07;       // infection rate
          double lambda = 0.; //2.68E-05;     // uninfected  birth rate
          double d = 0.; // 2.68E-05;   // uninfected death rate
          double psi = 2.82e-4; //.05;//1/365.; // rate to lose immunity (once a year is reasonable)
          double seasonality = 0.;
          double burnin = 0.; // percentage of events to discard

//        %$Pop 100 Locs 3 a 0.143 beta 0.0050 lambda=d 2.68E-5 time 60.0 migRate 0.0020$\\

// ----------------------------------------------------------------------------------------------------------------------------------------------------------------

		// Simulation parameters:
		int Ntraj = 1000;		// Number of trajectories
		int Nt = 100001;			// Number of timesteps
		int Nsamples = 101;	// Number of samples to record
		double T = 2000.;		// Length of time of simulation
		double alpha = 10;		// Critical reaction parameter

		// Model parameters:
        int s0 = 1000;
        int e0 = 0;
        int i0 = 1;
        int r0 = 0;

		double expose = 0.0;
		double infect = 8.14e-3/s0;
		double recover = 1.08e-3 + 2.82e-4;

		SEIRState x0 = new SEIRState(s0, e0, i0, r0, 0.0);

		// List to hold integrated trajectories:
		List<List<SEIRState>> trajectoryList = new ArrayList<List<SEIRState>>();

		// Initialize PRNG:
		Random random = new Random();
		// TODO: Need to set the seed for Poissonian RNG as well... How?

		// Create HybridTauleapSEIR instance:
		HybridTauleapSEIR hybridTauleapSEIR = new HybridTauleapSEIR(x0, expose, infect, recover, false, alpha, random);

        List<SEIRState> traj = hybridTauleapSEIR.genTrajectory(T, Nt, Nsamples, null);

        System.out.println("t\tS\tI\tR\n");


        for (SEIRState state : traj)
            System.out.println(state.time + "\t" + state.S + "\t" + state.I + "\t" + state.R);


		// Allocate and zero critical steps list:
		List<Integer> criticalTrajectories = new ArrayList<Integer>();
		for (int i=0; i<Nsamples; i++)
			criticalTrajectories.add(0);

		// Integrate trajectories:
		for (int i=0; i<Ntraj; i++) {
			hybridTauleapSEIR.setState(x0);
			trajectoryList.add(hybridTauleapSEIR.genTrajectory(T, Nt, Nsamples, criticalTrajectories));
		}



		// Calculate and print some moments:
//		List<SEIRStateDouble> means = SEIRStateMoments.getMeans(trajectoryList, Nsamples);
//		List<SEIRStateDouble> vars = SEIRStateMoments.getVariances(trajectoryList, means, Nsamples);
//		System.out.println("t S_mean S_var E_mean E_var I_mean I_var R_mean R_var cfrac");
//		for (int i=0; i<means.size(); i++) {
//			System.out.println(means.get(i).time + " "
//					+ means.get(i).S + " "
//					+ vars.get(i).S + " "
//					+ means.get(i).E + " "
//					+ vars.get(i).E + " "
//					+ means.get(i).I + " "
//					+ vars.get(i).I + " "
//					+ means.get(i).R + " "
//					+ vars.get(i).R + " "
//					+ (double)(criticalTrajectories.get(i))/Ntraj);
//		}


		// Done!
		System.exit(0);

	}
}
