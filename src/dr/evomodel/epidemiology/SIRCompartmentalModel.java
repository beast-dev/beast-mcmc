package dr.evomodel.epidemiology;

import dr.inference.model.Parameter;
import dr.math.MathUtils;     // for Randomizer.nextPoisson(), nextDouble()
import java.util.ArrayList;    // for CR and NCR lists
import java.util.List;
import cern.jet.random.Poisson;
import cern.jet.random.engine.RandomEngine;

public class SIRCompartmentalModel extends CompartmentalModel {

    Parameter transmissionRate;
    Parameter recoveryRate;
    Parameter samplingProportion;
    Parameter resusRate; // re-susceptibility rate
    Parameter numS;
    Parameter numI;
    Parameter numR;
    Parameter origin;
    Parameter tauStep;
    int numGridPoints;
    double cutOff;
    //Parameter epsilon;

    public SIRCompartmentalModel(
            Parameter transmissionRate,
            Parameter recoveryRate,
            Parameter samplingProportion,
            Parameter resusRate,
            Parameter numS,
            Parameter numI,
            Parameter numR,
            Parameter origin,
            Parameter tauStep,
            //Parameter epsilon,
            int numGridPoints,
            double cutOff) {

        super("SIRCompartmentalModel");

        this.transmissionRate = transmissionRate;
        addVariable(transmissionRate);
        this.recoveryRate = recoveryRate;
        addVariable(recoveryRate);
        this.samplingProportion = samplingProportion;
        addVariable(samplingProportion);
        this.resusRate = resusRate;
        addVariable(resusRate);
        this.numS = numS;
        addVariable(numS);
        this.numI = numI;
        addVariable(numI);
        this.numR = numR;
        addVariable(numR);
        this.origin = origin;
        addVariable(origin);
        this.tauStep = tauStep;
        addVariable(tauStep);
        //this.epsilon = epsilon;
        //addVariable(epsilon);

        this.numGridPoints = numGridPoints;
        this.cutOff = cutOff;
    }


    public void simulateTrajectory() {

        System.out.println("this is running");
        // print initial SIR
        System.out.println("initial numI: " + numI.getParameterValue(0));
        System.out.println("initial numR: " + numR.getParameterValue(0));
        System.out.println("initial numS: " + numS.getParameterValue(0));

        // set up time interval vector
        double T = origin.getParameterValue(0);
        double deltaTime = T / numGridPoints;
        double[] timeIntervals = new double[numGridPoints];
        for (int k = 0; k < numGridPoints; k++) {
            timeIntervals[k] = k * deltaTime;
        }

        // Initialize time and ensure index 0 is defined
        double time = 0.0;

        // already set initial values in numS/numI/numR and tauStep at index 0 via XML;
        // this reads them
        double S_i = numS.getParameterValue(0);
        double I_i = numI.getParameterValue(0);
        double R_i = numR.getParameterValue(0);

        // local trajectory storage on the fixed grid
        double[] S_traj = new double[numGridPoints];
        double[] I_traj = new double[numGridPoints];
        double[] R_traj = new double[numGridPoints];

        S_traj[0] = S_i;
        I_traj[0] = I_i;
        R_traj[0] = R_i;

        int nextRecordIndex = 1;

        // Simulate until we fill the grid (interval starts)
        while (time < T) {

            int K = numGridPoints;

            // 1) find current grid index i (0..K-1)
            int i = 0;
            for (int k = 0; k < K; k++) {
                if (timeIntervals[k] <= time) {
                    i = k;
                } else {
                    break;
                }
            }

            double beta = transmissionRate.getParameterValue(0);
            double gamma = recoveryRate.getParameterValue(0);
            double omega = resusRate.getParameterValue(0);


            // --------------------------------------------------------------------
            // TAU SELECTION
            // --------------------------------------------------------------------

            // default epsilon is 0.03
            double epsilon = 0.03;
            // a reaction is critical if there are fewer than criticalNumber left in that compartment
            // critical number usually between 2 and 20
            int criticalNumber = 5;

            // g_i's determined by using the highest order of reaction of species i
            double g_S = 2.0;
            double g_I = 2.0;
            double g_R = 1.0;

            double[] g = new double[]{g_S, g_I, g_R};

            // calculate epsilon using epsilon_i = epsilon / g_i
            double[] e = new double[]{
                    epsilon / g[0],
                    epsilon / g[1],
                    epsilon / g[2]
            };

            // step 1

            // v_ij
            double[] v_SI = new double[]{-1.0, 1.0, 0.0};
            double[] v_IR = new double[]{0.0, -1.0, 1.0};
            double[] v_RS = new double[]{1.0, 0.0, -1.0};

            // v matrix - now by row we have S, I, R
            double[][] vMatrix = new double[][]{
                    {-1.0, 0.0, 1.0},  // S
                    {1.0, -1.0, 0.0},  // I
                    {0.0, 1.0, -1.0}   // R
            };

            // Lj's: only one compartment is consumed in each reaction so Lj is only the one that is consumed
            double l_SI = Math.floor(S_i);
            double l_IR = Math.floor(I_i);
            double l_RS = Math.floor(R_i);

            boolean[] isCritical = new boolean[]{
                    l_SI < criticalNumber,
                    l_IR < criticalNumber,
                    l_RS < criticalNumber
            };

            List<Integer> cr = new ArrayList<>();
            List<Integer> ncr = new ArrayList<>();
            for (int j = 0; j < 3; j++) {
                if (isCritical[j]) cr.add(j);
                else ncr.add(j);
            }

            System.out.println("cr " + cr.size() + " ncr " + ncr.size());

            // step 2
            // tau selection only uses non-critical reactions

            // propensities (a_j)
            double a_SI = beta * S_i * I_i;
            double a_IR = gamma * I_i;
            double a_RS = omega * R_i;
            double a_0 = a_SI + a_IR + a_RS;

            double[] a = new double[]{a_SI, a_IR, a_RS};

            //check on a_0 <= 0
            if (a_0 <= 0.0) {
                while (nextRecordIndex < K) {
                    S_traj[nextRecordIndex] = S_i;
                    I_traj[nextRecordIndex] = I_i;
                    R_traj[nextRecordIndex] = R_i;
                    nextRecordIndex++;
                }
                break;
            }

            // calculate mu and sigma2
            double[] mu_hat = new double[3];
            double[] sigma2_hat = new double[3];

            for (int species = 0; species < 3; species++) {
                double mu = 0.0;
                double s2 = 0.0;

                for (int idx = 0; idx < ncr.size(); idx++) {
                    int reaction = ncr.get(idx);
                    double vij = vMatrix[species][reaction];
                    mu += vij * a[reaction];
                    s2 += (vij * vij) * a[reaction];
                }

                mu_hat[species] = mu;
                sigma2_hat[species] = s2;
            }

            // compute tau'

        /*
        // current counts
        double[] x = new double[]{S_i, I_i, R_i};

        // numerator
        double[] bound = new double[3];
        for (int species = 0; species < 3; species++) {
            bound[species] = Math.max(e[species] * x[species], 1.0);
        }

        double tau_prime = Double.POSITIVE_INFINITY;

        for (int species = 0; species < 3; species++) {

            double denom1 = Math.max(Math.abs(mu_hat[species]), 1e-16);
            double denom2 = Math.max(sigma2_hat[species], 1e-16);

            double tau1 = bound[species] / denom1;
            double tau2 = (bound[species] * bound[species]) / denom2;

            double candidate = Math.min(tau1, tau2);

            if (candidate < tau_prime) {
                tau_prime = candidate;
            }
        }
        */
            // compute tau'

            double tau_prime;

            //check this
            if (ncr.isEmpty()) {
                tau_prime = Double.POSITIVE_INFINITY;
            } else {

                // current counts
                double[] x = new double[]{S_i, I_i, R_i};

                // numerator
                double[] bound = new double[3];
                for (int species = 0; species < 3; species++) {
                    bound[species] = Math.max(e[species] * x[species], 1.0);
                }

                //initialize tau_prime
                tau_prime = Double.POSITIVE_INFINITY;

                for (int species = 0; species < 3; species++) {

                    // can't let denom be 0
                    double denom1 = Math.max(Math.abs(mu_hat[species]), 1e-16);
                    double denom2 = Math.max(sigma2_hat[species], 1e-16);

                    double tau1 = bound[species] / denom1;
                    double tau2 = (bound[species] * bound[species]) / denom2;

                    double candidate = Math.min(tau1, tau2);

                    if (candidate < tau_prime) {
                        tau_prime = candidate;
                    }
                }
            }

            // next grid boundary
            double nextBoundary = (nextRecordIndex < K) ? timeIntervals[nextRecordIndex] : T;

            // step 3

            boolean usedSSA = false;

            //typically use 10*1/a_0
            if (tau_prime < 10.0 / a_0) {

                //System.out.println("single step SAA started");

                //default is to run the single reaction SSA steps 100 times
                for (int k = 0; k < 100; k++) {
                    //calculate propensities
                    a_SI = beta * S_i * I_i;
                    a_IR = gamma * I_i;
                    a_RS = omega * R_i;
                    a_0 = a_SI + a_IR + a_RS;

                    if (a_0 == 0.0) {
                        break;
                    }

                    //find time to next reaction
                    double timeToRxn = -Math.log(MathUtils.nextDouble()) / a_0;

                    // if the next reaction is beyond the next grid boundary, just record current state there
                    if (time + timeToRxn >= nextBoundary) {
                        while (nextRecordIndex < K && timeIntervals[nextRecordIndex] <= Math.min(time + timeToRxn, T)) {
                            S_traj[nextRecordIndex] = S_i;
                            I_traj[nextRecordIndex] = I_i;
                            R_traj[nextRecordIndex] = R_i;
                            nextRecordIndex++;
                        }
                        time = Math.min(time + timeToRxn, T);
                        System.out.println("time: " + time);
                        break;
                    }

                    //choose reaction
                    //probability of each reaction occurring
                    double p_SI = a_SI / a_0;
                    double p_IR = a_IR / a_0;
                    double p_RS = a_RS / a_0;

                    double r = MathUtils.nextDouble();

                    int new_S;
                    int new_I;
                    int new_R;

                    if (r < p_SI) {
                        new_S = 0;
                        new_I = 1;
                        new_R = 0;
                    } else if (r < p_SI + p_IR) {
                        new_S = 0;
                        new_I = 0;
                        new_R = 1;
                    } else {
                        new_S = 1;
                        new_I = 0;
                        new_R = 0;
                    }

                    double S_new = S_i - new_I + new_S;
                    double I_new = I_i + new_I - new_R;
                    double R_new = R_i + new_R - new_S;

                    if (S_new < 0 || I_new < 1 || R_new < 0) {
                        break;
                    }

                    //update counts and time
                    //update time
                    time = time + timeToRxn;

                    System.out.println("time: " + time);

                    S_i = S_new;
                    I_i = I_new;
                    R_i = R_new;

                    // record current state at any grid points reached
                    while (nextRecordIndex < K && timeIntervals[nextRecordIndex] <= time + 1e-12) {
                        S_traj[nextRecordIndex] = S_i;
                        I_traj[nextRecordIndex] = I_i;
                        R_traj[nextRecordIndex] = R_i;
                        nextRecordIndex++;
                    }

                    //reset i
                    i = Math.max(0, nextRecordIndex - 1);

                    if (time >= T - 1e-12) {
                        break;
                    }

                    nextBoundary = (nextRecordIndex < K) ? timeIntervals[nextRecordIndex] : T;
                }

                System.out.println("initial numS: " + S_i + "numI: " + I_i + "numR: " + R_i);
                // if we do step 3, want to return to step 1 rather than proceed to step 4
                //System.out.println("SSA finished");
                usedSSA = true;
            }

            if (usedSSA) {
                continue;
            }

            // step 4

            // calculate a second option for tau
            double a0_critical = 0.0;
            for (int idx = 0; idx < cr.size(); idx++) {
                a0_critical += a[cr.get(idx)];
            }

            double tau_2prime;

            if (a0_critical > 0.0) {
                double u = MathUtils.nextDouble();
                tau_2prime = -Math.log(u) / a0_critical;
            } else {
                tau_2prime = Double.POSITIVE_INFINITY;
            }

            // step 5
            // determine tau
            double tau = Math.min(tau_prime, tau_2prime);

            // cap tau so we don't jump beyond the next grid boundary (keeps vectors consistent)
            // maybe this should go after user tau override so each boundary gets a count?
            tau = Math.min(tau, nextBoundary - time);
            tau = Math.min(tau, T - time);

            //this allows user to specify a tau step. default in xml should be 0
            if (tauStep.getParameterValue(0) > 0.0) {
                tau = tauStep.getParameterValue(0);
            }

            tau = Math.min(tau, nextBoundary - time);
            tau = Math.min(tau, T - time);

            System.out.println("tau:" + tau);

            //I don't remember what I was thinking with this. I'm pretty sure it stops if tau is too small and will take too long
            //if (tau <= cutOff) {
            //    System.out.println("breaking because tau <= cutOff");
            //    break;
            //}
            //does not work with this line of code

            boolean stepAccepted = false;
            double tau_try = tau;

            while (!stepAccepted) {

                if (tau_try <= 1e-16) {
                    break;
                }

                // --------------------------------------------------------------------
                // TAU LEAP
                // --------------------------------------------------------------------

                // assumes:
                // timeIntervals = interval starts, length K = numGridPoints
                // numS, numI, numR length K and already initialized at index 0
                // time is current time in [0, origin)
                // tau chosen elsewhere

                // 2) propensities
                a_SI = beta * S_i * I_i;
                a_IR = gamma * I_i;
                a_RS = omega * R_i;

                // 3) ODE derivatives for SAL (SIRS deterministic drift)
                double dSdt = -beta * S_i * I_i + omega * R_i;
                double dIdt = beta * S_i * I_i - gamma * I_i;
                double dRdt = gamma * I_i - omega * R_i;

                // 4) adot via chain rule
                double adot_SI = beta * (I_i * dSdt + S_i * dIdt);
                double adot_IR = gamma * dIdt;
                double adot_RS = omega * dRdt;

                // 5) SAL means
                double e_SI = a_SI * tau_try + 0.5 * adot_SI * tau_try * tau_try;
                double e_IR = a_IR * tau_try + 0.5 * adot_IR * tau_try * tau_try;
                double e_RS = a_RS * tau_try + 0.5 * adot_RS * tau_try * tau_try;

                e_SI = Math.max(e_SI, 0.0);
                e_IR = Math.max(e_IR, 0.0);
                e_RS = Math.max(e_RS, 0.0);

                // 6) sample events

            /*
            int new_I = Poisson.staticNextInt(e_SI);
            int new_R = Poisson.staticNextInt(e_IR);
            int new_S = Poisson.staticNextInt(e_RS);

            if (tau == tau_prime) {
                // for all critical reactions, no reaction occurs (step 5a)
                if (isCritical[0]) {
                    new_I = 0;
                }
                if (isCritical[1]) {
                    new_R = 0;
                }
                if (isCritical[2]) {
                    new_S = 0;
                }
            }

            // step 5b, one critical reaction occurs , other critical reactions are 0,
            if (tau == tau_2prime) {
                if (a0_critical > 0.0) {

                }
            }
            */

                // a more effective way of doing the above, this should be reusable with more reactions

                // 6) sample events

                double[] eVector = new double[]{e_SI, e_IR, e_RS};
                // reaction counts: [SI, IR, RS]
                int[] k_rxn = new int[3];

                if (tau_prime <= tau_2prime) {

                    // for all critical reactions, no reaction occurs (step 5a)
                    // for all noncritical reactions, generate new count as a Poisson sample.

                    for (int j = 0; j < 3; j++) {
                        if (isCritical[j]) {
                            k_rxn[j] = 0;
                        } else {
                            k_rxn[j] = Poisson.staticNextInt(eVector[j]);
                        }
                    }

                } else {

                    // step 5b, one critical reaction occurs, other critical reactions are 0,
                    // noncritical reactions are poisson sample
                    // the critical reaction that does occur is sampled same way as SSA

                    // first sample the noncritical reactions
                    for (int j = 0; j < 3; j++) {
                        if (isCritical[j]) {
                            k_rxn[j] = 0;
                        } else {
                            k_rxn[j] = Poisson.staticNextInt(eVector[j]);
                        }
                    }

                    // now choose jc among critical reactions only
                    if (a0_critical > 0.0) {
                        double[] a_rxn = new double[]{a_SI, a_IR, a_RS};

                        // draw uniform sample from (0, a0_critical)
                        double u = MathUtils.nextDouble() * a0_critical;
                        double cumulative = 0.0;
                        // jc is the index of the reaction that occurs, start at a value that it can't be
                        int jc = -1;

                        for (int idx = 0; idx < cr.size(); idx++) {
                            int j = cr.get(idx);
                            cumulative += a_rxn[j];
                            if (u < cumulative) {
                                jc = j;
                                break;
                            }
                        }

                        // exactly one critical reaction fires
                        if (jc >= 0) {
                            k_rxn[jc] = 1;
                        }
                    }
                }

                // number of new S, I, R's
                int new_I = k_rxn[0];
                int new_R = k_rxn[1];
                int new_S = k_rxn[2];

                if (new_R > I_i - 1) {
                    new_R = (int) Math.max(0, I_i - 1);
                }

                // 7) safety caps (debug-friendly; later replace with "shrink tau and redo")
                //new_I = Math.min(new_I, (int) Math.floor(S_i));
                //new_R = Math.min(new_R, (int) Math.floor(I_i + new_I));
                //new_S = Math.min(new_S, (int) Math.floor(R_i + new_R));

                // 8) update state after leap
                double S_new = S_i - new_I + new_S;
                double I_new = I_i + new_I - new_R;
                double R_new = R_i + new_R - new_S;

                // step 6 in paper - if any count is negative (or if less than 1 infected bc this means it won't run),
                // cut tow in half and go back to step 3
                if (S_new < 0 || I_new < 0 || R_new < 0) {
                    tau_try = tau_try / 2.0;
                    // continue should cause it to go back to step 3
                    continue;
                }

                // step 6 in paper, no count is negative so we continue
                stepAccepted = true;

                double time_new = time + tau_try;

                // 9) write to the correct grid index:
                // If you keep only K entries (interval starts), write to the interval containing time_new
                S_i = S_new;
                I_i = I_new;
                R_i = R_new;

                while (nextRecordIndex < K && timeIntervals[nextRecordIndex] <= time_new) {
                    S_traj[nextRecordIndex] = S_i;
                    I_traj[nextRecordIndex] = I_i;
                    R_traj[nextRecordIndex] = R_i;
                    nextRecordIndex++;
                }

                // print SIR
                System.out.println("numS: " + S_i + "numI: " + I_i + "numR: " + R_i);
                System.out.println("total: " + (S_i + I_i + R_i));

                time = time_new;

                // print time
                System.out.println("time: " + time);

                // stop if we’ve filled the last grid slot
                if (nextRecordIndex >= K) {
                    break;
                }
            }

            if (!stepAccepted) {
                break;
            }
        }

        // fill any remaining grid points with the last state
        while (nextRecordIndex < numGridPoints) {
            S_traj[nextRecordIndex] = S_i;
            I_traj[nextRecordIndex] = I_i;
            R_traj[nextRecordIndex] = R_i;
            nextRecordIndex++;
        }
    }
}

    // Code below is just for testing
        // for(int i = 0; i < numS.getDimension(); i++) {
        //     numS.setParameterValue(i, numS.getParameterValue(0) - i);
        // }

