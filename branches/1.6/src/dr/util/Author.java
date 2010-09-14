package dr.util;

/**
 * Container class of author infomation in a citation
 *
 * @author Alexei Drummond
 * @author Marc A. Suchard
 */

public class Author {

    String surname;
    String firstnames;

    public Author(String firstnames, String surname) {
        this.surname = surname;
        this.firstnames = firstnames;
    }

    public String getInitials() { // TODO Determine initials from first names
        return firstnames;
    }

    public String toString() {
        return surname + " " + getInitials();
    }
}
