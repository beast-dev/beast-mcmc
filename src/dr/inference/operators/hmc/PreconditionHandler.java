package dr.inference.operators.hmc;

public class PreconditionHandler {

    private final MassPreconditioner massPreconditioner;
    private final MassPreconditioningOptions options;
    private final MassPreconditionScheduler.Type schedulerType;

    public PreconditionHandler(MassPreconditioner massPreconditioner, MassPreconditioningOptions options, MassPreconditionScheduler.Type schedulerType) {
        this.massPreconditioner = massPreconditioner;
        this.options = options;
        this.schedulerType = schedulerType;
    }

    public MassPreconditioner getMassPreconditioner() {
        return massPreconditioner;
    }

    public MassPreconditioningOptions getOptions() {
        return options;
    }

    public MassPreconditionScheduler.Type getSchedulerType() {
        return schedulerType;
    }
}
