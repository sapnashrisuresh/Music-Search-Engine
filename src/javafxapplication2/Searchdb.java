/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package javafxapplication2;
/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */


import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.PreparedStatement;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.sql.Statement;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 *
 * @author nagashayan
 */
public class Searchdb {

    /**
     * @param args the command line arguments
     */
    //constructor to recieve mfcc values 
    public String mfcc="";
    private int sid;
    private String sname;
    private int gid;
    private String gname;
    private String album;
    private String artist,path;
    boolean flag=true;
    String arr[]=new String[5];  
    public Searchdb()
    {
         System.out.println("mfcc in search");
       /* for(int i=0;i<val.length;i++)
        {
            val[i]=Math.ceil(val[i]);
           mfcc=mfcc+val[i];
           System.out.println(mfcc); 
        }
          System.out.println(val.length);*/
        System.out.println("completed search"); 
    }
    
    public String[] search(int mm) throws SQLException
    {
       try
        {   
            Class.forName("com.mysql.jdbc.Driver");
             java.sql.Connection con=DriverManager.getConnection("jdbc:mysql://localhost:3306/pro_db?","root","");
             
             Statement stmt = null;
             String query ="select * " +
        "from pro_db.main1";
             stmt = con.createStatement();
        
            stmt.execute(query);
       

            ResultSet rs = stmt.executeQuery(query);
            while (rs.next()) {
             sid = rs.getInt("sid");
              sname = rs.getString("sname");
            //gid = rs.getInt("gid");
            gname = rs.getString("gname");
            album=rs.getString("album");
            artist=rs.getString("artist");
            path=rs.getString("path");
            System.out.println(sid + "\t" + sname +
                               "\t" + gid + "\t" + gname + "\t" +album + "\t" +artist + "\t" + path);
        
            //currently one max 15 columns are filled
        
            //for(int i=1;i<=5;i++)
          //{
            //  String m="";
              //m="m"+i;
              
              int match1=rs.getInt("m1");
              int match2=rs.getInt("m2");
              //int temp=Integer.parseInt(match)
            // int r=match.compareTo(mfcc);
            //boolean r1=match.equals(mm);
            if(match1==mm || match2==mm)
            {
                System.out.println("The song is"+sname);
                flag=false;
                arr[0]=sname;
            arr[1]=gname;
            arr[2]=album;
            arr[3]=artist;
            //arr[4]="C:\\Users\\nagashayan\\Documents\\NetBeansProjects\\JavaFXApplication2\\recorded.wav";
            arr[4]=path;
                System.out.println("The Song Name Is  "+arr[0]);
                
                
                break;
            }
            else if(mm<=(match2) && mm>=(match1))
            {
                System.out.println("The song is"+sname);
                flag=false;
                System.out.println("The song is"+sname);
                flag=false;
                arr[0]=sname;
            arr[1]=gname;
            arr[2]=album;
            arr[3]=artist;
           // arr[4]="C:\\Users\\nagashayan\\Documents\\NetBeansProjects\\JavaFXApplication2\\recorded.wav";
            arr[4]=path;
                System.out.println("The Song Name Is  "+arr[0]);
                
                break;
            }
            //else if()
          }
            //if(flag==false)
              // break;
        }   
            
             
           catch (ClassNotFoundException ex) {
            Logger.getLogger(Searchdb.class.getName()).log(Level.SEVERE, null, ex);
        }
        //}
       catch (SQLException ex) {
            Logger.getLogger(Searchdb.class.getName()).log(Level.SEVERE, null, ex);
        } 
        //}
        
        
        if(flag==true)
            {
                System.out.println("no match found");
              ///////////  String arr[] = null;
                arr[0]="Sorry, no match found!";
              //  return arr; //return "Sorry No Match Found!";
            }
            else
        {
            //String arr[] = null;
        /*    arr[0]=sname;
            arr[1]=gname;
            arr[2]=album;
            arr[3]=artist;
            arr[4]="C:\\Users\\nagashayan\\Documents\\NetBeansProjects\\JavaFXApplication2\\recorded.wav";
            
                System.out.println("The Song Name Is  "+arr[0]);*/
        //    return arr;
        }
       return arr;
    }//return sname+gname+album+artist;   
    }
    
    
    
    
    
    
    
    

      
