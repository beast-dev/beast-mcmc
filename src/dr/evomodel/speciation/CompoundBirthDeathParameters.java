package dr.evomodel.speciation;

import dr.evolution.tree.Tree;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;
import dr.inference.model.Bounds;
import dr.inference.model.VariableListener;

/**
 * Compound parameters for birth-death models supporting two modes:
 * 1. Compound-as-primary: R0, D, S are real parameters, λ, μ, ψ are proxies
 * 2. Raw-as-primary: λ, μ, ψ are real parameters, R0, D, S are proxies
 * 
 * Transformations:
 * R0 = λ/(μ + ψ)  - basic reproduction number
 * D = μ + ψ       - total removal rate
 * S = ψ/(μ + ψ)   - sampling proportion among removals
 * 
 * Inverse transformations:
 * λ = R0 * D
 * μ = D * (1 - S)  
 * ψ = D * S
 */
public class CompoundBirthDeathParameters {

    private final boolean compoundAsPrimary;
    
    // Primary parameters (the ones MCMC samples)
    private final Parameter primaryBirthRate;
    private final Parameter primaryDeathRate;
    private final Parameter primarySamplingRate;
    private final Parameter primaryR0;
    private final Parameter primaryD;
    private final Parameter primaryS;
    
    // Proxy parameters (computed from primary)
    public final Parameter birthRate;
    public final Parameter deathRate;
    public final Parameter samplingRate;
    public final Parameter R0;
    public final Parameter D;
    public final Parameter S;

    /**
     * Constructor for compound-as-primary mode (R0, D, S are the primary parameters)
     */
    public CompoundBirthDeathParameters(
            Parameter R0Parameter,
            Parameter DParameter,
            Parameter SParameter) {
            
        this.compoundAsPrimary = true;
        

        this.primaryR0 = R0Parameter;
        this.primaryD = DParameter;
        this.primaryS = SParameter;
        this.primaryBirthRate = null;
        this.primaryDeathRate = null;
        this.primarySamplingRate = null;
        

        this.birthRate = new BirthRateProxy("lambda", R0Parameter.getDimension());
        this.deathRate = new DeathRateProxy("mu", R0Parameter.getDimension());
        this.samplingRate = new SamplingRateProxy("psi", R0Parameter.getDimension());
        

        this.R0 = R0Parameter;
        this.D = DParameter;
        this.S = SParameter;
    }

    /**
     * Constructor for raw-as-primary mode (λ, μ, ψ are the primary parameters)
     */
    public CompoundBirthDeathParameters(
            Parameter birthRate,
            Parameter deathRate,
            Parameter samplingRate,
            boolean rawAsPrimary) {
            
        this.compoundAsPrimary = false;
        
        // Store primary raw parameters
        this.primaryBirthRate = birthRate;
        this.primaryDeathRate = deathRate;
        this.primarySamplingRate = samplingRate;
        this.primaryR0 = null;
        this.primaryD = null;
        this.primaryS = null;
        

        this.birthRate = birthRate;
        this.deathRate = deathRate;
        this.samplingRate = samplingRate;
        

        this.R0 = new R0Proxy("R0", birthRate.getDimension());
        this.D = new DProxy("D", birthRate.getDimension());
        this.S = new SProxy("S", birthRate.getDimension());
    }

    /**
     * Proxy for λ = R0 * D (when compound parameters are primary)
     */
    private class BirthRateProxy extends Parameter.Proxy implements VariableListener {

        public BirthRateProxy(String name, int dimension) {
            super(name, dimension);
            primaryR0.addVariableListener(this);
            primaryD.addVariableListener(this);
        }

        @Override
        public double getParameterValue(int i) {
            return primaryR0.getParameterValue(i) * primaryD.getParameterValue(i);
        }

        @Override
        public void setParameterValue(int i, double value) {

        }

        @Override
        public void setParameterValueQuietly(int i, double value) {

        }

        @Override
        public void setParameterValueNotifyChangedAll(int i, double value) {

        }

        @Override
        public void variableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
            fireParameterChangedEvent();
        }
    }

    /**
     * Proxy for μ = D * (1 - S) (when compound parameters are primary)
     */
    private class DeathRateProxy extends Parameter.Proxy implements VariableListener {

        public DeathRateProxy(String name, int dimension) {
            super(name, dimension);
            primaryD.addVariableListener(this);
            primaryS.addVariableListener(this);
        }

        @Override
        public double getParameterValue(int i) {
            return primaryD.getParameterValue(i) * (1.0 - primaryS.getParameterValue(i));
        }

        @Override
        public void setParameterValue(int i, double value) {

        }

        @Override
        public void setParameterValueQuietly(int i, double value) {

        }

        @Override
        public void setParameterValueNotifyChangedAll(int i, double value) {

        }

        @Override
        public void variableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
            fireParameterChangedEvent();
        }
    }

    /**
     * Proxy for ψ = D * S (when compound parameters are primary)
     */
    private class SamplingRateProxy extends Parameter.Proxy implements VariableListener {

        public SamplingRateProxy(String name, int dimension) {
            super(name, dimension);
            primaryD.addVariableListener(this);
            primaryS.addVariableListener(this);
        }

        @Override
        public double getParameterValue(int i) {
            return primaryD.getParameterValue(i) * primaryS.getParameterValue(i);
        }

        @Override
        public void setParameterValue(int i, double value) {

        }

        @Override
        public void setParameterValueQuietly(int i, double value) {

        }

        @Override
        public void setParameterValueNotifyChangedAll(int i, double value) {

        }

        @Override
        public void variableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
            fireParameterChangedEvent();
        }
    }

    /**
     * Proxy for R0 = λ/(μ + ψ) (when raw parameters are primary)
     */
    private class R0Proxy extends Parameter.Proxy implements VariableListener {

        public R0Proxy(String name, int dimension) {
            super(name, dimension);
            primaryBirthRate.addVariableListener(this);
            primaryDeathRate.addVariableListener(this);
            primarySamplingRate.addVariableListener(this);
        }

        @Override
        public double getParameterValue(int i) {
            double totalRemoval = primaryDeathRate.getParameterValue(i) + primarySamplingRate.getParameterValue(i);
            return primaryBirthRate.getParameterValue(i) / totalRemoval;
        }

        @Override
        public void setParameterValue(int i, double value) {

        }

        @Override
        public void setParameterValueQuietly(int i, double value) {

        }

        @Override
        public void setParameterValueNotifyChangedAll(int i, double value) {

        }

        @Override
        public void variableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
            fireParameterChangedEvent();
        }
    }

    /**
     * Proxy for D = μ + ψ (when raw parameters are primary)
     */
    private class DProxy extends Parameter.Proxy implements VariableListener {

        public DProxy(String name, int dimension) {
            super(name, dimension);
            primaryDeathRate.addVariableListener(this);
            primarySamplingRate.addVariableListener(this);
        }

        @Override
        public double getParameterValue(int i) {
            return primaryDeathRate.getParameterValue(i) + primarySamplingRate.getParameterValue(i);
        }

        @Override
        public void setParameterValue(int i, double value) {

        }

        @Override
        public void setParameterValueQuietly(int i, double value) {

        }

        @Override
        public void setParameterValueNotifyChangedAll(int i, double value) {

        }

        @Override
        public void variableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
            fireParameterChangedEvent();
        }
    }

    /**
     * Proxy for S = ψ/(μ + ψ) (when raw parameters are primary)
     */
    private class SProxy extends Parameter.Proxy implements VariableListener {

        public SProxy(String name, int dimension) {
            super(name, dimension);
            primaryDeathRate.addVariableListener(this);
            primarySamplingRate.addVariableListener(this);
        }

        @Override
        public double getParameterValue(int i) {
            double totalRemoval = primaryDeathRate.getParameterValue(i) + primarySamplingRate.getParameterValue(i);
            return primarySamplingRate.getParameterValue(i) / totalRemoval;
        }

        @Override
        public void setParameterValue(int i, double value) {

        }

        @Override
        public void setParameterValueQuietly(int i, double value) {

        }

        @Override
        public void setParameterValueNotifyChangedAll(int i, double value) {

        }

        @Override
        public void variableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
            fireParameterChangedEvent();
        }
    }

    // Getter methods
    public Parameter getBirthRate() {
        return birthRate;
    }

    public Parameter getDeathRate() {
        return deathRate;
    }

    public Parameter getSamplingRate() {
        return samplingRate;
    }

    public Parameter getR0Parameter() {
        return R0;
    }

    public Parameter getDParameter() {
        return D;
    }

    public Parameter getSParameter() {
        return S;
    }

    public boolean isCompoundAsPrimary() {
        return compoundAsPrimary;
    }

    public boolean isCompoundParameter(Variable variable) {
        return variable == R0 || variable == D || variable == S;
    }

    public boolean isRawParameter(Variable variable) {
        return variable == birthRate || variable == deathRate || variable == samplingRate;
    }
} 