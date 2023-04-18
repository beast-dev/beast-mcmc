package dr.inference.hmc;


public abstract class IncrementTransformType {
    IncrementTransformType(String transformType) {
        this.transformType = transformType;
    }

    public String getTransformType() {
        return transformType;
    }

    private String transformType;
    public abstract double getDerivativeOfInverseTransform(double x);
    public abstract double[] parameterFromIncrements(double[] x);
    public abstract double[] incrementsFromParameter(double[] x);

    public static IncrementTransformType factory(String match, double upper, double lower){
        if (match.equalsIgnoreCase("log")) {
            return new LogTransform(match);
        } else if (match.equalsIgnoreCase("logit")) {
            return new LogitTransform(match, upper, lower);
        }
        return null;
    }
}

class LogTransform extends IncrementTransformType {
    LogTransform(String transformType) {
        super(transformType);
    }

    public double getDerivativeOfInverseTransform(double x) {
        return Math.exp(x);
    }

    public double[] parameterFromIncrements(double[] delta) {
        double[] fx = new double[delta.length];
        fx[0] = delta[0];
        for (int i = 1; i < delta.length; i++) {
            fx[i] = fx[i-1] + delta[i];
        }
        for (int i = 0; i < delta.length; i++) {
            fx[i] = Math.exp(fx[i]);
        }
        return fx;
    }

    public double[] incrementsFromParameter(double[] x) {
        double[] increments = new double[x.length];
        increments[0] = Math.log(x[0]);
        for (int i = 1; i < x.length; i++) {
            increments[i] = Math.log(x[i]/x[i - 1]);
        }
        return increments;
    }
}

class LogitTransform extends IncrementTransformType {
    LogitTransform(String transformType, double upper, double lower) {
        super(transformType);
        this.upper = upper;
        this.lower = lower;
    }
    private double upper = Double.POSITIVE_INFINITY;
    private double lower = Double.NEGATIVE_INFINITY;

    private double getScaledLogit(double x) {
        double u = (x - lower / (upper - lower));
        return Math.log(u / (1 - u));
    }

    private double getScaledSigmoid(double y) {
        double inverse = 1 / (1 + Math.exp(-y));
        return lower + (upper-lower) * inverse;
    }
    public double getDerivativeOfInverseTransform(double y) {
        double inverse = 1 / (1 + Math.exp(-y));
        return (upper-lower) * inverse * (1 - inverse);
    }

    public double[] parameterFromIncrements(double[] delta) {
        double[] fx = new double[delta.length];
        fx[0] = delta[0];
        for (int i = 1; i < delta.length; i++) {
            fx[i] = fx[i-1] + delta[i];
        }
        for (int i = 0; i < delta.length; i++) {
            fx[i] = getScaledSigmoid(fx[i]);
        }
        return fx;
    }

    public double[] incrementsFromParameter(double[] x) {
        double[] increments = new double[x.length];
        increments[0] = getScaledLogit(x[0]);
        for (int i = 1; i < x.length; i++) {
            increments[i] = getScaledLogit(x[i]) - getScaledLogit(x[i-1]);
        }
        return increments;
    }

    public double getUpper() {return upper;}

    public double getLower() {return lower;}

}
