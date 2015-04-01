/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package javafxapplication2;

import java.io.File;
import java.io.IOException;
import java.net.URL;
import java.util.ResourceBundle;
import java.util.Scanner;
import javafx.beans.value.ChangeListener;
import javafx.beans.value.ObservableValue;
import javafx.concurrent.Task;
import javafx.fxml.FXML;
import javafx.fxml.Initializable;
import javafx.scene.control.Button;
import javafx.scene.control.Label;
import javafx.scene.control.ProgressIndicator;
import javax.sound.sampled.AudioFileFormat;
import javax.sound.sampled.AudioFormat;
import javax.sound.sampled.AudioInputStream;
import javax.sound.sampled.AudioSystem;
import javax.sound.sampled.Clip;
import javax.sound.sampled.DataLine;
import javax.sound.sampled.LineUnavailableException;
import javax.sound.sampled.SourceDataLine;
import javax.sound.sampled.TargetDataLine;
import javax.sound.sampled.UnsupportedAudioFileException;

//=============================================//
//Inner class to play back the data from the
// audio file.
/**
 *
 * @author nagashayan
 */
public class UIController implements Initializable {
    
    AudioFormat audioFormat;
      TargetDataLine targetDataLine;
       AudioFormat audioFormat1;
  AudioInputStream audioInputStream;
  SourceDataLine sourceDataLine;
  String path;
      
    @FXML
    private Button stop;
    @FXML
    private Button record;
    @FXML
    private Button stop1;
    @FXML
    private Button record1;
     @FXML
    private Button stop2;
    @FXML
    private Label label1;
    @FXML
    final ProgressIndicator progress = new ProgressIndicator(0);
    int xrand;
    Task searchworker;
    Task searchworker2;
    @FXML
    private void stop2()
    {
        stop2.setDisable(true);
        //targetDataLine.stop();
          //targetDataLine.close();
   //       new playsong().start();
        //label.setText("Music Search Engine");
          //progress.setProgress(0);
          searchworker2 = createWorker1();
          progress.progressProperty().unbind();
          progress.progressProperty().bind(searchworker2.progressProperty()); 
       
      searchworker2.messageProperty().addListener(new ChangeListener<String>() {
                    public void changed(ObservableValue<? extends String> observable, String oldValue, String newValue) {
                        System.out.println("buddy the song is"+newValue);
                        label1.setText(newValue);
                    }
                });

                new Thread(searchworker2).start();
         
         
        
    }
    public Task createWorker1() {
      return new Task() {
          
            @Override
            protected Object call() throws Exception {
                     //label.setText("Music Search Engine");
          //Matlab obj=new Matlab();
          //obj.main2();
          updateProgress(1, 10);
          //progress.setProgress(0.1);
          //MFCC IMPLEMENTATION GOES HEREEEEEEEEEEEEE.....
           
             
             
          
          int nnumberofFilters = 24; 
          int nlifteringCoefficient = 22; 
          boolean oisLifteringEnabled = true; 
          boolean oisZeroThCepstralCoefficientCalculated = true; 
          int nnumberOfMFCCParameters = 12; //without considering 0-th 
          double dsamplingFrequency = 44100.0; 
          int nFFTLength = 524288; 
          updateProgress(2, 10);
          if (oisZeroThCepstralCoefficientCalculated) { 
            //take in account the zero-th MFCC 
            nnumberOfMFCCParameters = nnumberOfMFCCParameters + 1; 
          } 
//          else { 
//            nnumberOfMFCCParameters = nnumberOfMFCCParameters; 
//          } 
          updateProgress(3, 10);
          MFCC mfcc = new MFCC(nnumberOfMFCCParameters, 
                               dsamplingFrequency, 
                               nnumberofFilters, 
                               nFFTLength, 
                               oisLifteringEnabled, 
                               nlifteringCoefficient, 
                               oisZeroThCepstralCoefficientCalculated); 
          isDebugMode = true;
 
          debug(mfcc.toString()); 
          updateProgress(4, 10);
          double[] x = new double[16777216]; 
          Double[] y=new Double[524288];
          double num; 
          try{     	  
        	  int i;
                  String rand="C:\\Users\\Home\\Documents\\NetBeansProjects\\JavaFXApplication2\\clips\\"+xrand+".txt";
          Scanner fileScan = new Scanner (new File(rand));
          for (i=0; i<5242880 ; i++){
        	  
        	  if (fileScan.hasNext())
        	  {
        	  num = fileScan.nextDouble();
        	  x[i] = num;
        	  
        	  }
          }
          
          }
         
          catch (Exception e)
          {
              e.printStackTrace();
              System.out.println("No such file exists.");
          }
          Double[] dparameters = null;
          //int m=0;
          int mm = 0;
          double big=0,small;
          int temp;
          String mfc="";
          for(int k=0;k<524288;k+=524288)
          {
        	 // System.out.println("loop 1");
        	  int j=0;
          	for (int i=k; i<k+524288; i++){
        	  y[j]=x[i];
        	  j++;
          }
            dparameters = mfcc.getParameters(y); 
          //debug("MFCC parameters:"); 
            small=dparameters[0];
          for (int i = 0; i < (dparameters.length); i++) { 
        	//  debug(" " + dparameters[i]); 
              System.out.println(dparameters[i]);
              temp=(int) Math.round(dparameters[i]);
              double temp1=Math.round(dparameters[i]);
              mfc=mfc+temp1;
              big=Math.max(big,temp);
              small=Math.min(small, temp);
          }
           System.out.println("maximum is:"+big);
          System.out.println("minimum is:"+small);
          mm=     (int) (big-small);
          big=0;
          
          small=0;
          System.out.println("max-min is:"+mm);
          //System.out.println("mean is:"+a1[2]/2);
            //m++;   
          } 	 
         
            
        
        updateProgress(5, 10);  
          System.out.println(dparameters.length);
          //just for ease of running use above code
         
            System.out.println("completed mfcc");
            
           //Searchdb search=new Searchdb(dparameters);
           updateProgress(6, 10);
           
            RandSearchdb search=new RandSearchdb(mfc);
            //search.search(mm);
            String[] arr=search.search();
          String sname=arr[0];
          sname="Song name: "+sname+"\n"+"Genre:"+arr[1]+"\n"+"Album:"+arr[2]+"\n"+"Artist:"+arr[3];
        // String sname=arr[0];
         System.out.println(sname);
         //path=arr[4];
         //System.out.println(path);
           updateProgress(7, 10);
          //System.out.println("The Song Name is "+sname);
         
          updateProgress(8, 10);
          updateProgress(9, 10);
          updateMessage(sname);
                    updateProgress(10, 10);
                   //label1.setText(sname);
               System.out.println(mfc);    
               stop2.setDisable(false);
         return true;

            
                }
              };
          
    }
    @FXML
    private void record1()
    {
        
    }
    @FXML
    private void randplay()
    {
        xrand=(int)Math.round((Math.random()*10));
        String rand="C:\\Users\\Home\\Documents\\NetBeansProjects\\JavaFXApplication2\\clips\\"+xrand+".wav";
        try{
      File soundFile =new File(rand);
    
        //int x=Math.random();
      //  File soundFile =new File(path);
        audioInputStream = AudioSystem.
                  getAudioInputStream(soundFile);
      audioFormat1 = audioInputStream.getFormat();
      System.out.println(audioFormat1);

      DataLine.Info dataLineInfo =
                          new DataLine.Info(
                            SourceDataLine.class,
                                    audioFormat1);

      sourceDataLine =(SourceDataLine)AudioSystem.getLine(dataLineInfo);

     
      new RandPlayThread().start();
    }catch (Exception e) {
      e.printStackTrace();
      System.exit(0);
    }

    }
    class RandPlayThread extends Thread{
  byte tempBuffer[] = new byte[10000];
    

  public void run(){
      //boolean check=stop1.isDisabled();
      //System.out.println(check);
      
      try
      {
      sourceDataLine.open(audioFormat1);
      sourceDataLine.start();

      int cnt;
      
      while((cnt = audioInputStream.read(
           tempBuffer,0,tempBuffer.length)) != -1 && (stop2.isDisabled()==false)){
        if(cnt > 0){
          //Write data to the internal buffer of
          // the data line where it will be
          // delivered to the speaker.
          sourceDataLine.write(
                             tempBuffer, 0, cnt);
        }//end if
      }//end while
      //Block and wait for internal buffer of the
      // data line to empty.
      sourceDataLine.drain();
      sourceDataLine.close();

      //Prepare to playback another file
      
    }catch (Exception e) {
      e.printStackTrace();
      System.exit(0);
    }//end catch
  }//end run
}//end inner class PlayThread
//===================================//
    @FXML
    private void stop1() throws IOException
    {
        stop1.setDisable(true);
    }
    

    
    @FXML
    private void play1() throws IOException
    {
      
        System.out.println("am in play event");
       
    try{
      //File soundFile =new File("C:\\Users\\nagashayan\\Documents\\NetBeansProjects\\JavaFXApplication2\\recorded2.wav");
      File soundFile =new File(path);
        audioInputStream = AudioSystem.
                  getAudioInputStream(soundFile);
      audioFormat1 = audioInputStream.getFormat();
      System.out.println(audioFormat1);

      DataLine.Info dataLineInfo =
                          new DataLine.Info(
                            SourceDataLine.class,
                                    audioFormat1);

      sourceDataLine =(SourceDataLine)AudioSystem.getLine(dataLineInfo);

     
      new PlayThread().start();
    }catch (Exception e) {
      e.printStackTrace();
      System.exit(0);
    }//end catch
  }//end play()

//===================================//


    class PlayThread extends Thread{
  byte tempBuffer[] = new byte[10000];
    

  public void run(){
      //boolean check=stop1.isDisabled();
      //System.out.println(check);
      
      try
      {
      sourceDataLine.open(audioFormat1);
      sourceDataLine.start();

      int cnt;
      
      while((cnt = audioInputStream.read(
           tempBuffer,0,tempBuffer.length)) != -1 && (stop1.isDisabled()==false)){
        if(cnt > 0){
          //Write data to the internal buffer of
          // the data line where it will be
          // delivered to the speaker.
          sourceDataLine.write(
                             tempBuffer, 0, cnt);
        }//end if
      }//end while
      //Block and wait for internal buffer of the
      // data line to empty.
      sourceDataLine.drain();
      sourceDataLine.close();

      //Prepare to playback another file
      
    }catch (Exception e) {
      e.printStackTrace();
      System.exit(0);
    }//end catch
  }//end run
}//end inner class PlayThread
//===================================//


    
    @FXML
    private void recordButtonAction() {
                
        //System.out.println("You clicked record!");
        label1.setText("Record for 15-20 seconds...");
        record.setDisable(true);
        stop.setDisable(false);
         captureAudio();
        //label.setText("Music Search Engine");
         
               
    }
    
    @FXML
     private void stopButtonAction() {
        
       // System.out.println("You clicked stop!");
       
      record.setDisable(false);
      stop.setDisable(true);
         targetDataLine.stop();
          targetDataLine.close();
   //       new playsong().start();
        //label.setText("Music Search Engine");
          //progress.setProgress(0);
          searchworker = createWorker();
          progress.progressProperty().unbind();
          progress.progressProperty().bind(searchworker.progressProperty()); 
       
       searchworker.messageProperty().addListener(new ChangeListener<String>() {
                    public void changed(ObservableValue<? extends String> observable, String oldValue, String newValue) {
                        System.out.println("buddy the song is"+newValue);
                        label1.setText(newValue);
                    }
                });

                new Thread(searchworker).start();
         
                         
    }
     public static boolean isDebugMode;
       public static void debug(String str){
        	if(isDebugMode){
        		System.out.println("[" + MFCC.class.getSimpleName()+ "]"+str);
        	}
        }
    //create a Task search worker
    
    public Task createWorker() {
      return new Task() {
          
            @Override
            protected Object call() throws Exception {
                     //label.setText("Music Search Engine");
          Matlab obj=new Matlab();
          obj.main2();
          updateProgress(1, 10);
          //progress.setProgress(0.1);
          //MFCC IMPLEMENTATION GOES HEREEEEEEEEEEEEE.....
           
             
             
          
          int nnumberofFilters = 24; 
          int nlifteringCoefficient = 22; 
          boolean oisLifteringEnabled = true; 
          boolean oisZeroThCepstralCoefficientCalculated = true; 
          int nnumberOfMFCCParameters = 12; //without considering 0-th 
          double dsamplingFrequency = 44100.0; 
          int nFFTLength = 524288; 
          updateProgress(2, 10);
          if (oisZeroThCepstralCoefficientCalculated) { 
            //take in account the zero-th MFCC 
            nnumberOfMFCCParameters = nnumberOfMFCCParameters + 1; 
          } 
//          else { 
//            nnumberOfMFCCParameters = nnumberOfMFCCParameters; 
//          } 
          updateProgress(3, 10);
          MFCC mfcc = new MFCC(nnumberOfMFCCParameters, 
                               dsamplingFrequency, 
                               nnumberofFilters, 
                               nFFTLength, 
                               oisLifteringEnabled, 
                               nlifteringCoefficient, 
                               oisZeroThCepstralCoefficientCalculated); 
          isDebugMode = true;
 
          debug(mfcc.toString()); 
          updateProgress(4, 10);
          double[] x = new double[16777216]; 
          Double[] y=new Double[524288];
          double num; 
          try{     	  
        	  int i;
          Scanner fileScan = new Scanner (new File("amprecorded.txt"));
          for (i=0; i<5242880 ; i++){
        	  
        	  if (fileScan.hasNext())
        	  {
        	  num = fileScan.nextDouble();
        	  x[i] = num;
        	  
        	  }
          }
          
          }
         
          catch (Exception e)
          {
              e.printStackTrace();
              System.out.println("No such file exists.");
          }
          Double[] dparameters = null;
          //int m=0;
          int mm = 0;
          double big=0,small;
          int temp;
          for(int k=0;k<524288;k+=524288)
          {
        	 // System.out.println("loop 1");
        	  int j=0;
          	for (int i=k; i<k+524288; i++){
        	  y[j]=x[i];
        	  j++;
          }
            dparameters = mfcc.getParameters(y); 
          //debug("MFCC parameters:"); 
            small=dparameters[0];
          for (int i = 0; i < (dparameters.length-1); i++) { 
        	//  debug(" " + dparameters[i]); 
              System.out.println(dparameters[i]);
              temp=(int) Math.round(dparameters[i]);
              
              big=Math.max(big,temp);
              small=Math.min(small, temp);
          }
           System.out.println("maximum is:"+big);
          System.out.println("minimum is:"+small);
          mm=     (int) (big-small);
          big=0;
          
          small=0;
          System.out.println("max-min is:"+mm);
          //System.out.println("mean is:"+a1[2]/2);
            //m++;   
          } 	 
         
            
        
        updateProgress(5, 10);  
          System.out.println(dparameters.length);
          //just for ease of running use above code
         
            System.out.println("completed mfcc");
            
           //Searchdb search=new Searchdb(dparameters);
           updateProgress(6, 10);
           
            Searchdb search=new Searchdb();
            //search.search(mm);
            
          String[] arr=search.search(mm);
          
            
          String sname=arr[0];
          sname="Song name: "+sname+"\n"+"Genre:"+arr[1]+"\n"+"Album:"+arr[2]+"\n"+"Artist:"+arr[3];
         //String sname=arr[0];
         System.out.println(sname);
         path=arr[4];
         System.out.println(path);
           updateProgress(7, 10);
          //System.out.println("The Song Name is "+sname);
         
          updateProgress(8, 10);
          updateProgress(9, 10);
          updateMessage(sname);
                    updateProgress(10, 10);
                   //label1.setText(sname);
                   
         return true;

            
                }
              };
          
    }  
          
     
    
            

            

          
    
    
    @FXML
    private void onMouseEnterR(){
        record.setScaleX(1.3);
         record.setScaleY(1.3);
             
    }
    
    @FXML
    private void onMouseExitR(){
        record.setScaleX(1);
         record.setScaleY(1);
     
    }
    
    @FXML
    private void onMouseEnterS(){
        stop.setScaleX(1.3);
         stop.setScaleY(1.3);
     
    }
    
    @FXML
    private void onMouseExitS(){
        stop.setScaleX(1);
         stop.setScaleY(1);
     
    }
    
    @Override
    public void initialize(URL url, ResourceBundle rb) {
        // TODO
      
        
        
        
    }  
    
    private void captureAudio(){
    try{
     
      audioFormat = getAudioFormat();
      DataLine.Info dataLineInfo =
                          new DataLine.Info(
                            TargetDataLine.class,
                            audioFormat);
      targetDataLine = (TargetDataLine)
               AudioSystem.getLine(dataLineInfo);

     
      new CaptureThread().start();
    }catch (Exception e) {
      e.printStackTrace();
      System.exit(0);
    }
    
    
    
  }

    private AudioFormat getAudioFormat(){
    float sampleRate = 44100.0F;
   
    int sampleSizeInBits = 16;
   
    int channels = 2;
   
    boolean signed = true;
   
    boolean bigEndian = false;

    return new AudioFormat(sampleRate,
                           sampleSizeInBits,
                           channels,
                           signed,
                           bigEndian);
  }

    class CaptureThread extends Thread{
  public void run(){
    AudioFileFormat.Type fileType = null;
    File audioFile = null;
      fileType = AudioFileFormat.Type.WAVE;
      audioFile = new File("recorded.wav");
   

    try{
      targetDataLine.open(audioFormat);
      targetDataLine.start();
      AudioSystem.write(
            new AudioInputStream(targetDataLine),
            fileType,
            audioFile);
      //JTextArea area=new JTextArea();
      //area.setText(targetDataLine);
    }catch (Exception e){
      e.printStackTrace();
    }

  }
}
   /* class playsong extends Thread
    {
        public void run()
        {
            AudioInputStream audioIn = null;
            try {
                URL url = this.getClass().getClassLoader().getResource("recorded2.wav");
                audioIn = AudioSystem.getAudioInputStream(url);
                // Get a sound clip resource.
                Clip clip = AudioSystem.getClip();
                // Open audio clip and load samples from the audio input stream.
                clip.open(audioIn);
                clip.start();
            } catch (UnsupportedAudioFileException ex) {
                Logger.getLogger(UIController.class.getName()).log(Level.SEVERE, null, ex);
            } catch (IOException ex) {
                Logger.getLogger(UIController.class.getName()).log(Level.SEVERE, null, ex);
            } catch (LineUnavailableException ex) {
                Logger.getLogger(UIController.class.getName()).log(Level.SEVERE, null, ex);
            } finally {
                try {
                    audioIn.close();
                } catch (IOException ex) {
                    Logger.getLogger(UIController.class.getName()).log(Level.SEVERE, null, ex);
                }
            }

        }
    }*/
}
            
