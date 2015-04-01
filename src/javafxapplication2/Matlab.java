package javafxapplication2;

//package javafxapplication2;

import java.io.*;

public class Matlab extends JavaFXApplication2
{
  
  public static String runScript(File f1) {
	  String output = " ", error = " ";
	  try {
		 // String path="C:\\Users\\Home\\Desktop\\matlab\\trial.m";
              
	  String commandToRun = "C:\\Program Files (x86)\\MATLAB\\R2009a\\bin\\matlab.exe /C START /MIN -nojvm -nodisplay -nosplash -nodesktop -wait -r run('C:\\Users\\Home\\Documents\\NetBeansProjects\\JavaFXApplication2\\trial.m');exit;";
	  //"C:\<a long path here>\matlab.exe" -nodisplay -nosplash -nodesktop -r "run('C:\<a long path here>\mfile.m');"
	  System.out.println(commandToRun);
	  Process p = Runtime.getRuntime().exec(commandToRun);

	  String s;

	  BufferedReader stdInput = new BufferedReader(new
	  InputStreamReader(p.getInputStream()));

	  BufferedReader stdError = new BufferedReader(new
	  InputStreamReader(p.getErrorStream()));

	  // read the output from the command
	 // System.out.println("\nHere is the standard output of the command:\n");
	  while ((s = stdInput.readLine()) != null) {
	  output += s + "\n";
	  System.out.println(s);
	  }

	  // read any errors from the attempted command
	  //System.out.println("\nHere is the standard error of the command (if any):\n");
	  while ((s = stdError.readLine()) != null) {
	  error += s + "\n";
	  System.out.println(s);
	  }

	  } catch (Exception e) {
	  //System.out.println("exception happened � here�s what I know:");
	  e.printStackTrace();
	  System.exit(-1);
	  }
	  return output + error;
	  }
  
  public int main2()
  {
	  String path="C:\\Users\\Home\\Desktop\\matlab\\trial.m";
	  File f1=new File(path);
	  runScript(f1);
	  return 1;
	
  }
}