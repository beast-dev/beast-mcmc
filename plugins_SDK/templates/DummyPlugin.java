package @PLUGIN_FULL_PACKAGE@;
import java.util.*;

import @PLUGIN_FULL_PACKAGE@.*;
import dr.app.plugin.*;
import dr.xml.XMLObjectParser;
public class @PLUGIN_CLASS@ implements Plugin {

	public Set<XMLObjectParser> getParsers() {
		Set<XMLObjectParser> parsers = new HashSet<XMLObjectParser>();
		@PLUGIN_CLASS@Parser pluginParser = new @PLUGIN_CLASS@Parser();
		parsers.add(pluginParser);
		return parsers;
	}
	
}
