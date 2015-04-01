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
public class RandSearchdb {

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
    private String artist;
    boolean flag=true;
    String arr[]=new String[5];  
    public RandSearchdb(String val)
    {
         System.out.println("mfcc in search");
        //for(int i=0;i<val.length;i++)
        //{
            //val[i]=Math.ceil(val[i]);
          // mfcc=mfcc+val[i];
            mfcc=val;
           System.out.println(mfcc); 
        //}
          //System.out.println(val.length);
        System.out.println("completed search"); 
    }
    
    public String[] search()
    {
       try
        {   
            Class.forName("com.mysql.jdbc.Driver");
             java.sql.Connection con=DriverManager.getConnection("jdbc:mysql://localhost:3306/pro_db?","root","");
             
             Statement stmt = null;
             String query ="select * " +
        "from pro_db.main2";
             stmt = con.createStatement();
        
            stmt.execute(query);
       

            ResultSet rs = stmt.executeQuery(query);
            while (rs.next()) {
             sid = rs.getInt("sid");
              sname = rs.getString("sname");
           // gid = rs.getInt("gid");
            gname = rs.getString("gname");
            album=rs.getString("album");
            artist=rs.getString("artist");
            System.out.println(sid + "\t" + sname +
                                "\t" + gname + "\t" +album + "\t" +artist);
        
            //currently one max 15 columns are filled
        
            for(int i=1;i<=5;i++)
          {
              String m="";
              m="m"+i;
              
              //int match=rs.getInt(m);
              String match=rs.getString(m);
              //int temp=Integer.parseInt(match)
            // int r=match.compareTo(mfcc);
            //boolean r1=match.equals(mm);
            //String match=rs.getString(m);
          //String mfcc="0256.088.0-269.0-29.0227.0198.0214.0-335.0-405.0180.0225.0";
           
             
            //int r=match.compareTo(mfcc);
            boolean r1=match.equals(mfcc);
            if(r1==true)
            {
                System.out.println("The song is"+sname);
                flag=false;
                arr[0]=sname;
            arr[1]=gname;
            arr[2]=album;
            arr[3]=artist;
                
                break;
            }
          }
            if(flag==false)
               break;
            
            
            }
           
        } catch (SQLException ex) {
            Logger.getLogger(Searchdb.class.getName()).log(Level.SEVERE, null, ex);
            
        } catch (ClassNotFoundException ex) {
            Logger.getLogger(Searchdb.class.getName()).log(Level.SEVERE, null, ex);
        }
        if(flag==true)
            {
                System.out.println("no match found");
                arr[0]="Sorry, no match found!";
            arr[1]="";
            arr[2]="";
            arr[3]="";
                 
            }
            else
        {
                //return ("The Song Name Is  "+sname);
       
            arr[0]=sname;
            arr[1]=gname;
            arr[2]=album;
            arr[3]=artist;
        } 
        return arr;
       //return sname+gname+album+artist;   
    }
}
    
    
    
    
    
    
    
  