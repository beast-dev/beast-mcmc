package dr.util;

import java.io.*;

/**
 * @author joseph
 *         Date: 15/04/2008
 */
public class FileHelpers {

    public static final String FILE_NAME = "fileName";

    /**
     * @param file the file to read the line numbers from
     * @return Number of lines in file
     * @throws IOException low level file error
     */
    public static int numberOfLines(File file) throws IOException {
        RandomAccessFile randFile = new RandomAccessFile(file, "r");
        long lastRec = randFile.length();
        randFile.close();
        FileReader fileRead = new FileReader(file);
        LineNumberReader lineRead = new LineNumberReader(fileRead);
        lineRead.skip(lastRec);
        int count = lineRead.getLineNumber() - 1;
        fileRead.close();
        lineRead.close();
        return count;
    }

    /**
     * Resolve file from name.
     * <p/>
     * Keep A fully qualified (i.e. absolute path) as is. A name starting with a "./" is
     * relative to the master directory (set by FileHelpers.setMasterDir).
     * Any other name is stripped of any directory
     * component and placed in the "user.dir" directory.
     *
     * @param fileName an absolute or relative file name
     * @return a File object resolved from provided file name
     */
    public static File getFile(String fileName) {
        final boolean localFile = fileName.startsWith("./");
        final boolean relative = masterDirectory != null && localFile;
        if (localFile) {
            fileName = fileName.substring(2);
        }

        final File file = new File(fileName);
        final String name = file.getName();
        String parent = file.getParent();

        if (!file.isAbsolute()) {
            String p;
            if (relative) {
                p = masterDirectory.getAbsolutePath();
            } else {
                p = System.getProperty("user.dir");
            }
            if (parent != null && parent.length() > 0) {
                parent = p + '/' + parent;
            } else {
                parent = p;
            }
        }
        return new File(parent, name);
    }


    // directory where beast xml file resides
    private static File masterDirectory = null;

    public static void setMasterDir(File fileName) {
        masterDirectory = fileName;
    }
}
