package dr.app.plugin;

import dr.xml.XMLObjectParser;
import java.util.*;

public interface Plugin {
    Set<XMLObjectParser> getParsers();
}
