package dr.app.bss;

import java.awt.Frame;
import java.awt.Image;
import java.awt.Toolkit;
import java.net.URL;

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
	
	public Image CreateImage(String path) {
		URL imgURL = this.getClass().getResource(path);
		Toolkit kit = Toolkit.getDefaultToolkit();
		Image img = kit.createImage(imgURL);

		if (img != null) {
			return img;
		} else {
			System.err.println("Couldn't find file: " + path + "\n");
			return null;
		}

	}// END: CreateImage
	
}// END: class
