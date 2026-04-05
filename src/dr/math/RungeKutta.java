package dr.math;

/**
 * Fixed-step classical fourth-order Runge-Kutta (RK4) integrator for systems of ODEs.
 *
 * Solves dy/dt = f(t, y) using pre-allocated scratch arrays to avoid per-step allocation.
 *
 * @author Frederik M. Andersen
 */
public class RungeKutta {

    private final double[] k1;
    private final double[] k2;
    private final double[] k3;
    private final double[] k4;
    private final double[] tmp;

    @FunctionalInterface
    public interface RhsFunction {
        void evaluate(double t, double[] y, double[] dydt);
    }

    /**
     * @param systemSize dimension of the ODE system
     */
    public RungeKutta(int systemSize) {
        this.k1 = new double[systemSize];
        this.k2 = new double[systemSize];
        this.k3 = new double[systemSize];
        this.k4 = new double[systemSize];
        this.tmp = new double[systemSize];
    }

    /**
     * Single RK4 step: advance y(t) -> y(t+h).
     *
     * @param t    current time
     * @param h    step size
     * @param y    input state vector
     * @param yOut output state vector (may alias y)
     * @param n    number of equations to integrate
     * @param rhs  right-hand side function
     */
    public void step(double t, double h, double[] y, double[] yOut, int n, RhsFunction rhs) {
        final double h2 = 0.5 * h;
        final double h6 = h / 6.0;
        final double h3 = h / 3.0;

        rhs.evaluate(t, y, k1);
        for (int i = 0; i < n; i++) tmp[i] = y[i] + h2 * k1[i];

        rhs.evaluate(t + h2, tmp, k2);
        for (int i = 0; i < n; i++) tmp[i] = y[i] + h2 * k2[i];

        rhs.evaluate(t + h2, tmp, k3);
        for (int i = 0; i < n; i++) tmp[i] = y[i] + h * k3[i];

        rhs.evaluate(t + h, tmp, k4);
        for (int i = 0; i < n; i++) {
            yOut[i] = y[i] + h6 * (k1[i] + k4[i]) + h3 * (k2[i] + k3[i]);
        }
    }
}
