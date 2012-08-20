package dr.util;

import java.util.List;

/**
 * Interface for associating a list of citations with an object
 *
 * @author Marc A. Suchard
 */

public interface Citable {

    /**
     * @return a list of citations associated with this object
     */
    List<Citation> getCitations();

    public class Utils {

        public static String getCitationString(Citable citable, String prepend, String postpend) {
            List<Citation> citations = citable.getCitations();
            if (citations == null || citations.size() == 0) {
                return null;
            }
            StringBuilder builder = new StringBuilder();
            for (Citation citation : citations) {
                builder.append(prepend);
                builder.append(citation.toString());
                builder.append(postpend);
            }
            return builder.toString();
        }

        public static String getCitationString(Citable citable) {
            return getCitationString(citable, DEFAULT_PREPEND, DEFAULT_POSTPEND);
        }

        public static final String DEFAULT_PREPEND = "\t\t";
        public static final String DEFAULT_POSTPEND = "\n";
    }
}
