package dr.app.bss;

import java.lang.Thread.UncaughtExceptionHandler;
import java.net.URL;

import javax.swing.ImageIcon;
import javax.swing.JOptionPane;
import javax.swing.SwingUtilities;

public class ExceptionHandler implements UncaughtExceptionHandler {

	ImageIcon errorIcon = CreateImageIcon("icons/error.png");

	public void uncaughtException(final Thread t, final Throwable e) {

		if (SwingUtilities.isEventDispatchThread()) {
			showException(t, e);
		} else {
			SwingUtilities.invokeLater(new Runnable() {
				public void run() {
					showException(t, e);
				}
			});
		}
	}

	private void showException(Thread t, Throwable e) {
		String msg = String.format("Unexpected problem on thread %s: %s", t
				.getName(), e.getMessage());

		logException(t, e);

		JOptionPane.showMessageDialog(Utils.getActiveFrame(), msg, "Error",
				JOptionPane.ERROR_MESSAGE, errorIcon);
	}

	private void logException(Thread t, Throwable e) {
		// TODO: start a thread that logs it
		e.printStackTrace();

	}

	private ImageIcon CreateImageIcon(String path) {
		URL imgURL = this.getClass().getResource(path);
		if (imgURL != null) {
			return new ImageIcon(imgURL);
		} else {
			System.err.println("Couldn't find file: \n" + path + "\n");
			return null;
		}
	}

}
