package dr.evomodel.epidemiology;

public class SEIRState {

    public int S;
    public int E;
    public int I;
    public int R;
    public double time;

    public SEIRState(int S, int E, int I, int R, double time) {

        this.S = S;
        this.E = E;
        this.I = I;
        this.R = R;
        this.time = time;
    }

    public SEIRState copy() {
    	return new SEIRState(S, E, I, R, time);
    }
}
