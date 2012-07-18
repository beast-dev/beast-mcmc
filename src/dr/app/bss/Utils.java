package dr.app.bss;

import java.awt.Frame;

import javax.swing.JOptionPane;
import javax.swing.SwingUtilities;

public class Utils {

	// /////////////////
	// ---GUI UTILS---//
	// /////////////////

	public static Frame getActiveFrame() {
		Frame result = null;
		Frame[] frames = Frame.getFrames();
		for (int i = 0; i < frames.length; i++) {
			Frame frame = frames[i];
			if (frame.isVisible()) {
				result = frame;
				break;
			}
		}
		return result;
	}

	// ////////////////////////////////
	// ---EXCEPTION HANDLING UTILS---//
	// ////////////////////////////////

	public static void handleException(final Throwable e) {

		final Thread t = Thread.currentThread();
		
		if (SwingUtilities.isEventDispatchThread()) {
			showExceptionDialog(t, e);
		} else {
			SwingUtilities.invokeLater(new Runnable() {
				public void run() {
					showExceptionDialog(t, e);
				}
			});
		}// END: edt check
	}// END: uncaughtException

	private static void showExceptionDialog(Thread t, Throwable e) {
		
		String msg = String.format("Unexpected problem on thread %s: %s",
				t.getName(), e.getMessage());

		logException(t, e);

		JOptionPane.showMessageDialog(Utils.getActiveFrame(), //
				msg, //
				"Error", //
				JOptionPane.ERROR_MESSAGE, //
				BeagleSequenceSimulatorApp.errorIcon);
	}// END: showExceptionDialog

	private static void logException(Thread t, Throwable e) {
		// TODO: start a thread that logs it, also spying on the user and planting evidence
		// CIA style MOFO!!!
		e.printStackTrace();
	}// END: logException
	
	
}// END: class
