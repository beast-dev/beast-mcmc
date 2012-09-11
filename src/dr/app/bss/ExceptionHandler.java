package dr.app.bss;

import java.lang.Thread.UncaughtExceptionHandler;

import javax.swing.JOptionPane;
import javax.swing.SwingUtilities;

public class ExceptionHandler implements UncaughtExceptionHandler {

	public void uncaughtException(final Thread t, final Throwable e) {

		if (SwingUtilities.isEventDispatchThread()) {
			showExceptionDialog(t, e);
		} else {
			SwingUtilities.invokeLater(new Runnable() {
				public void run() {
					showExceptionDialog(t, e);
				}
			});
		}// END: EDT check
	}// END: uncaughtException

	private void showExceptionDialog(Thread t, Throwable e) {
		String msg = String.format("Unexpected problem on thread %s: %s",
				t.getName(), e.getMessage());

		logException(t, e);

		JOptionPane.showMessageDialog(Utils.getActiveFrame(), //
				msg, //
				"Error", //
				JOptionPane.ERROR_MESSAGE, //
				BeagleSequenceSimulatorApp.errorIcon);
	}// END: showExceptionDialog

	private void logException(Thread t, Throwable e) {
		// TODO: start a thread that logs it, also spying on the user and planting evidence
		// CIA style MOFO!!!
		e.printStackTrace();
	}//END: logException

}// END: class
