package dr.inference.operators.hmc;

public class PreconditionHandler {

    private final MassPreconditioner type;
    private final MassPreconditioningOptions options;
    private final MassPreconditionScheduler scheduler;

    public PreconditionHandler(MassPreconditioner type, MassPreconditioningOptions options, MassPreconditionScheduler scheduler) {
        this.type = type;
        this.options = options;
        this.scheduler = scheduler;
    }

    public MassPreconditioner getType() {
        return type;
    }

    public MassPreconditioningOptions getOptions() {
        return options;
    }

    public MassPreconditionScheduler getScheduler() {
        return scheduler;
    }
}
