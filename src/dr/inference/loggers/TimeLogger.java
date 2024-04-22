package dr.inference.loggers;

import dr.util.Timer;

public class TimeLogger implements Loggable {

    private final Timer timer;
    private Boolean hasStarted = false;

    public TimeLogger() {
        this.timer = new Timer();
    }


    @Override
    public LogColumn[] getColumns() {
        return new LogColumn[]{
                new LogColumn.Abstract("secondsElapsed") {

                    @Override
                    protected String getFormattedValue() {
                        if (!hasStarted) {
                            hasStarted = true;
                            timer.start();
                        }

                        return Double.toString(timer.toSeconds());
                    }
                }
        };
    }
}
