package dr.util;

/**
 * @author Alexei Drummond
 */
public class Citation {

    Author[] authors;
    String title;
    int year;
    String journal;
    int volume;
    int startpage;
    int endpage;

    public Citation() {
    }

    public String toString() {

        StringBuilder builder = new StringBuilder();
        builder.append("<html>");
        builder.append(authors[0].toString());
        for (int i = 1; i < authors.length; i++) {
            builder.append(", ");
            builder.append(authors[i].toString());
        }
        builder.append(" (").append(year).append(") ");
        builder.append(title).append(". ");
        builder.append("<i>").append(journal).append("</i>");
        builder.append(" <b>").append(volume).append("</b>:");
        builder.append(startpage).append("-").append(endpage);
        builder.append("</html>");

        return builder.toString();
    }

    class Author {

        String surname;
        String firstnames;

        public String getInitials() {
            return firstnames;
        }

        public String toString() {

            return surname + " " + getInitials();

        }

    }

}
