package dr.app.plugin;

import java.io.File;
import java.io.FileFilter;
import java.net.URL;
import java.net.URLClassLoader;
import java.util.ArrayList;
import java.util.List;
import java.util.logging.Logger;

public class PluginLoader {

	  
	   static final String PLUGIN_FOLDER = "plugins";
	   static final File PLUGIN_FILE = new File(PLUGIN_FOLDER);

//	    private static Set<String> loadedPlugins = new HashSet<String>();

	   public static List<String> getAvailablePlugins(){

	       List<String> plugins = new ArrayList<String> ();
	       File[] classFolderFiles = PLUGIN_FILE.listFiles(new FileFilter() {
	           public boolean accept(File pathname) {
	               String name = pathname.getName();
	               if(!pathname.isDirectory() || name.endsWith("CVS") || name.endsWith(".classes"))
	                   return false;
	               File[] directoryContents = pathname.listFiles(new FileFilter() {
	                   public boolean accept(File pathname) {
	                       String name = pathname.getName();
	                       return name.endsWith(".jar");
	                   }
	               });
	               return directoryContents.length != 0;
	           }
	       });

	       if (classFolderFiles != null) {
	           for (File folder : classFolderFiles) {
	               plugins.add(folder.getName());
	           }
	       }

	       File[] pluginJarFiles = PLUGIN_FILE.listFiles(new FileFilter() {
	           public boolean accept(File pathname) {
	               return !pathname.isDirectory() && pathname.getAbsolutePath().endsWith(".jar");
	           }
	       });

	       if (pluginJarFiles != null) {
	           for (File jarFile : pluginJarFiles) {
	               String name = jarFile.getName();
	               name =name.substring(0, name.length()- 4);
	               if(! plugins.contains(name))
	                   plugins.add(name);
	           }
	       }

	       return plugins;
	   }

	  public static Plugin loadPlugin(final String pluginName/*, boolean pluginEnabled*/) {   //the class loader must still be assigned if the plugin isnt enabled so
	                                                                                     //documents from that plugin can still be displayed.
          final String loggerName = "dr.app.plugin";
          Logger.getLogger(loggerName).info("loading plugin " + pluginName);
	      String fullname = PLUGIN_FOLDER + "/" + pluginName;
	      File file = new File(fullname);

	      try {
	          URL[] urls;
	          if (!file.exists()) {
	        	  Logger.getLogger(loggerName).info("loading jar file");
	              fullname = PLUGIN_FOLDER + "/" + pluginName+ ".jar";
	              file = new File(fullname);
	              urls = new URL[]{file.toURL()};
	          }
	          else {
	              File classFiles = new File(fullname + "/" + "classes");
	              final boolean classesExist = classFiles.exists();

	              File[] files = file.listFiles(new FileFilter() {
	                  public boolean accept(File pathname) {
	                      String name = pathname.getName();
	                      if(!name.endsWith(".jar")) return false;
	                      name = name.substring(0, name.length()- 4);
	                      return !(name.equals(pluginName) && classesExist);
	                  }
	              });
	              if(files == null) files = new File[0];

	              int length = files.length+1;
	              if(classesExist) length ++;

	              urls = new URL[length];
	              int count = 0;
	              if (classesExist) {
	                  urls [ count ++ ] =classFiles.toURL();
	              }
	              urls [ count ++ ] = file.toURL();
	              Logger.getLogger(loggerName).info("adding " + file + " to class path");
	              for (File jarFile : files) {
	                  urls[count++] = jarFile.toURL();
	              }
	          }
	          URLClassLoader classLoader;
	          classLoader = new URLClassLoader(urls);

	          for (URL url : classLoader.getURLs()) {
	        	  Logger.getLogger(loggerName).info("URL from loader: " + url.toString() + "\n");
	          }
	          
	          Class myClass = classLoader.loadClass(pluginName);
	          
	          Object plugin = myClass.newInstance();
	          if (!(plugin instanceof Plugin)){
	              throw new Exception("Class should be " + Plugin.class.getName());
	          }
	          return (Plugin)plugin;
	          
	          
	      } catch (Exception e) {
	    	  Logger.getLogger(loggerName).severe(e.getMessage());
	      }
	      return null;
	  }
}
