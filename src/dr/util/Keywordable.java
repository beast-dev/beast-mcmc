package dr.util;

import java.util.List;

/**
 * @author Andrew Rambaut
 * @version $Id$
 */
public interface Keywordable {

    void addKeyword(String keyword);

    List<String> getKeywords();

}
