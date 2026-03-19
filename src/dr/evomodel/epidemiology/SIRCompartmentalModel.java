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
    Parameter taustep;
    int numGridPoints;
    double cutOff;

    public SIRCompartmentalModel(
            Parameter transmissionRate,
            Parameter recoveryRate,
            Parameter samplingProportion,
            Parameter resusRate,
            Parameter numS,
            Parameter numI,
            Parameter numR,
            Parameter origin,
            Parameter taustep,
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
        this.taustep = taustep;
        addVariable(taustep);

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

        // already set initial values in numS/numI/numR and taustep at index 0 via XML;
        // this reads them
        numS.setParameterValue(0, numS.getParameterValue(0));
        numI.setParameterValue(0, numI.getParameterValue(0));
        numR.setParameterValue(0, numR.getParameterValue(0));

        // Simulate until we fill the grid (interval starts)
        while (time < T - 1e-12) {

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

            double S_i = numS.getParameterValue(i);
            double I_i = numI.getParameterValue(i);
            double R_i = numR.getParameterValue(i);

            double beta = transmissionRate.getParameterValue(0);
            double gamma = recoveryRate.getParameterValue(0);
            double omega = resusRate.getParameterValue(0);

            // --------------------------------------------------------------------
            // TAU SELECTION
            // --------------------------------------------------------------------

            // default epsilon
            double epsilon = 0.03;
            int criticalNumber = 10;

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
            double[][] V = new double[][]{
                    {-1.0, 0.0, 1.0},  // S
                    {1.0, -1.0, 0.0},  // I
                    {0.0, 1.0, -1.0}   // R
            };

            // Lj's: only one compartment is consumed in each reaction so Lj is only the one that is consumed
            double L_SI = Math.floor(S_i);
            double L_IR = Math.floor(I_i);
            double L_RS = Math.floor(R_i);

            boolean[] isCritical = new boolean[]{
                    L_SI < criticalNumber,
                    L_IR < criticalNumber,
                    L_RS < criticalNumber
            };

            List<Integer> CR = new ArrayList<>();
            List<Integer> NCR = new ArrayList<>();
            for (int j = 0; j < 3; j++) {
                if (isCritical[j]) CR.add(j);
                else NCR.add(j);
            }

            // step 2
            // tau selection only uses non-critical reactions

            // propensities (a_j)
            double a_SI = beta * S_i * I_i;
            double a_IR = gamma * I_i;
            double a_RS = omega * R_i;

            double[] a = new double[]{a_SI, a_IR, a_RS};

            // calculate mu and sigma2
            double[] mu_hat = new double[3];
            double[] sigma2_hat = new double[3];

            for (int species = 0; species < 3; species++) {
                double mu = 0.0;
                double s2 = 0.0;

                for (int idx = 0; idx < NCR.size(); idx++) {
                    int reaction = NCR.get(idx);
                    double vij = V[species][reaction];
                    mu += vij * a[reaction];
                    s2 += (vij * vij) * a[reaction];
                }

                mu_hat[species] = mu;
                sigma2_hat[species] = s2;
            }

            // compute tau'

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

            // step 3

            // step 4

            // calculate a second option for tau
            double a0_critical = 0.0;
            for (int idx = 0; idx < CR.size(); idx++) {
                a0_critical += a[CR.get(idx)];
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
            double nextBoundary = (i < K - 1) ? timeIntervals[i + 1] : T;
            tau = Math.min(tau, nextBoundary - time);
            tau = Math.min(tau, T - time);

            //this allows user to specify a tau step. default in xml should be 0
            if (taustep.getParameterValue(0) > 0.0) {
                tau = taustep.getParameterValue(0);
            }

            System.out.println("tau:" + tau);

            //I don't remember what I was thinking with this. I'm pretty sure it stops if tau is too small and will take too long
            //if (tau <= cutOff) {
            //    System.out.println("breaking because tau <= cutOff");
            //    break;
            //}

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
            double E_SI = a_SI * tau + 0.5 * adot_SI * tau * tau;
            double E_IR = a_IR * tau + 0.5 * adot_IR * tau * tau;
            double E_RS = a_RS * tau + 0.5 * adot_RS * tau * tau;

            E_SI = Math.max(E_SI, 0.0);
            E_IR = Math.max(E_IR, 0.0);
            E_RS = Math.max(E_RS, 0.0);

            // 6) sample events
            int new_I = Poisson.staticNextInt(E_SI);
            int new_R = Poisson.staticNextInt(E_IR);
            int new_S = Poisson.staticNextInt(E_RS);

            // 7) safety caps (debug-friendly; later replace with "shrink tau and redo")
            new_I = Math.min(new_I, (int) Math.floor(S_i));
            new_R = Math.min(new_R, (int) Math.floor(I_i + new_I));
            new_S = Math.min(new_S, (int) Math.floor(R_i + new_R));

            // 8) update state after leap
            double S_new = S_i - new_I + new_S;
            double I_new = I_i + new_I - new_R;
            double R_new = R_i + new_R - new_S;

            double time_new = time + tau;

            // 9) write to the correct grid index:
            // If you keep only K entries (interval starts), write to the interval containing time_new
            double tForIndex = Math.min(time_new, T - 1e-12);

            int j = 0;
            for (int k = 0; k < K; k++) {
                if (timeIntervals[k] <= tForIndex) {
                    j = k;
                } else {
                    break;
                }
            }

            numS.setParameterValue(j, S_new);
            numI.setParameterValue(j, I_new);
            numR.setParameterValue(j, R_new);

            // print SIR
            System.out.println("numS: " + S_new + "numI: " + I_new + "numR: " + R_new);

            time = time_new;

            // print time
            System.out.println("time: " + time);

            // stop if we’ve filled the last grid slot
            if (j >= K - 1) {
                break;
            }
        }
    }
}

    // Code below is just for testing
        // for(int i = 0; i < numS.getDimension(); i++) {
        //     numS.setParameterValue(i, numS.getParameterValue(0) - i);
        // }

