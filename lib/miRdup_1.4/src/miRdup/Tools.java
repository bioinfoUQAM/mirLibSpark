/*
 *  miRdup v1.0
 *  Computational prediction of the localization of microRNAs within their pre-miRNA
 *  
 *  Copyright (C) 2013  Mickael Leclercq
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 * 
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the 
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

/*
 * tools
 */
package miRdup;

import java.io.*;
import java.lang.reflect.Method;
import java.net.MalformedURLException;
import java.net.URL;
import java.net.URLClassLoader;
import java.net.URLConnection;
import java.util.Enumeration;
import java.util.zip.ZipEntry;
import java.util.zip.ZipFile;


/**
 *
 * @author Mickael Leclercq
 */
public class Tools {
    
     /**
     * Get HTML code from an url
     * @param link
     * @return code
     */
    public static String getCode (String link) {
        if (Main.debug) System.out.println("Getting source code of "+link);
        String sourceCode = null;
        try {
            URL url = new URL(link);
            URLConnection uc = url.openConnection();
            InputStream in = uc.getInputStream();
            int c = in.read();
            StringBuilder build = new StringBuilder();
            while (c != -1) {
                build.append((char) c);
                c = in.read();
            }
            sourceCode = build.toString();
        } catch (MalformedURLException e) {
            e.printStackTrace();
        } catch (IOException e) {
            e.printStackTrace();
        }
        return sourceCode;
    }
    
    public static void saveToFile(String s, String outfile){
        if (Main.debug) System.out.println("Saving "+outfile);
        try {
            PrintWriter pw = new PrintWriter(new FileWriter(outfile,false));
            pw.println(s);
            pw.close();
        } catch (Exception e) {
            saveToFile(s, outfile);
        }
    }
    
    /**
     * Download a file from an url
     * @param weblink 
     */
    public static void downloadFile(String weblink){
        if (Main.debug) System.out.println("Downloading "+weblink);
        try {
            URL url = new URL(weblink);
            File outfile= new File(weblink.substring(weblink.lastIndexOf("/")+1));
            BufferedInputStream bis = new BufferedInputStream(url.openStream());
            FileOutputStream fos = new FileOutputStream(outfile);
            BufferedOutputStream bos = new BufferedOutputStream(fos, 1024);
            byte[] data = new byte[1024];
            int i=0;
            while((i=bis.read(data,0,1024))>=0){
                bos.write(data,0,i);
            }
            bos.close();
            bis.close();
        } catch (Exception e) {
            downloadFile(weblink);
        }
    }
    
    /**
     * Unzip a file
     * @param zippedFilename 
     */
    public static void unzipFile(String zippedFilename){
        if (Main.debug) System.out.println("Unzipping "+zippedFilename);
        Enumeration entries;
        ZipFile zipFile;
        try {
            zipFile = new ZipFile(zippedFilename);
            entries = zipFile.entries();
            while(entries.hasMoreElements()) {
                ZipEntry entry = (ZipEntry)entries.nextElement();
                if(entry.isDirectory()) {
                    // Assume directories are stored parents first then children.
                    System.err.println("Extracting directory: " + entry.getName());
                    // This is not robust, just for demonstration purposes.
                    (new File(entry.getName())).mkdir();
                    continue;
                }
                try {
                    InputStream is = zipFile.getInputStream(entry);
                    OutputStream os = new BufferedOutputStream(new FileOutputStream(entry.getName()));
                    byte[] buffer = new byte[1024];
                    int len;
                    while((len = is.read(buffer)) >= 0)
                    os.write(buffer, 0, len);
                    is.close();
                    os.close();
                } catch (Exception e) {
                    e.printStackTrace();
                }
            }
            zipFile.close();
        } catch (IOException ioe) {
            System.err.println("Unhandled exception:");
            ioe.printStackTrace();
            return;
        }
    }
    
    /**
     * Delete a file
     * @param file 
     */
    public static void deleteFile(String file){
        if (Main.debug) System.out.println("Deleting "+file);
        File f = new File(file);
        f.delete();
    }
    
    /**
     * Rename a file
     * @param oldName
     * @param newName 
     */
    public static void renameFile(String oldName, String newName){
        if (Main.debug) System.out.println("Renaming "+oldName+" in "+newName);
        File old= new File(oldName);
        File newf= new File(newName);
        old.renameTo(newf);
    }
    
    /**
     * Execute linux command with bash
     * @param cmd
     * @return 
     */
    public static String executeLinuxCommand (String cmd){
        if (Main.debug) System.out.println("Executing command "+cmd);
        try {
            ProcessBuilder pb = new ProcessBuilder("bash", "-c", cmd);
            pb.redirectErrorStream(true); // use this to capture messages sent to stderr
            Process shell = pb.start();
            InputStream is = shell.getInputStream(); // this captures the output from the command
            int shellExitStatus = shell.waitFor(); 
            
            StringBuilder response = new StringBuilder();
              int value = 0;
              boolean active = true;
              while (active) {
                    value = is.read();
                    
                    if (value == -1) {
                        active=false;
                    } else {
                        response.append((char) value);
                        continue;
                    }
                }
                return response.toString();
            
        } catch (Exception e) {
            return "";
        }
    }
       
    /**
     * Execute linux command with bash 
     * @param cmd
     * @return 
     */
    public static void executeLinuxCommandToFile (String cmd){
        if (Main.debug) System.out.println("Executing command "+cmd);
        try {
            ProcessBuilder pb = new ProcessBuilder("bash", "-c", cmd);
            pb.redirectErrorStream(true); // use this to capture messages sent to stderr
            Process shell = pb.start();
            InputStream is = shell.getInputStream(); // this captures the output from the command
            int shellExitStatus = shell.waitFor(); 
            
            StringBuilder response = new StringBuilder();
              int value = 0;
              boolean active = true;
              while (active) {
                    value = is.read();
                    if (value == -1) {
                        active = false;
                        //throw new IOException("End of Stream");
                    } else  {
                        response.append((char) value);
                        continue;
                    } 
                }
              if (!response.toString().trim().isEmpty()) {
                  System.err.println("RNAfold says: "+response);
              }
            
        } catch (Exception e) {
            e.printStackTrace();
        }
    }
     
     
     /**
     * Execute windows command with cmd
     * @param cmd
     */
    public static String executeWindowsCommand2(String cmd) {
        //String path=getcmdpath().trim();
        Runtime runtime = Runtime.getRuntime();
        String output="";
        try{
            //Process ps=pb.start();
            String c = "cmd /c "+cmd;
            Process ps=runtime.exec(c);
            if (Main.debug) System.out.println("Executing command "+c);
            BufferedReader bri=new BufferedReader(new InputStreamReader(ps.getInputStream()));
            BufferedReader bre = new BufferedReader(new InputStreamReader(ps.getErrorStream()));
            String line="";
            while ((line=bri.readLine())!=null) {
               //System.out.println(line);
               output=output+line+"\n";
            }
            while ((line = bre.readLine()) != null) {
               System.out.println(line);
            }
            bre.close();
            ps.waitFor();
        }
        catch(Exception e){
            e.printStackTrace();
        }
        return output;
    }
    
         /**
         * Execute windows command with cmd
         * @param cmd
         */
    public static String executeWindowsCommand(String cmd) {

        //System.out.println(cmd);
        Runtime runtime = Runtime.getRuntime();
        String output="";
        try{
            //Process ps=pb.start();
            Process ps=runtime.exec("cmd /C "+cmd);
            
            BufferedReader bri=new BufferedReader(new InputStreamReader(ps.getInputStream()));
            String line="";
            while ((line=bri.readLine())!=null) {
               //System.out.println(line);
               output=output+line+"\n";
           }
           ps.waitFor();
        }
        catch(Exception e){
            e.printStackTrace();
        }
        return output;
    }
     /**
     * Execute windows command with cmd
     * @param cmd
     */
    public static void executeWindowsCommandToFile(String cmd) {
        //String path=getcmdpath().trim();
        Runtime runtime = Runtime.getRuntime();
        try{
            //Process ps=pb.start();
            String c = "cmd /c "+cmd;
            Process ps=runtime.exec(c);
            if (Main.debug) System.out.println("Executing command "+c);
            BufferedReader bri=new BufferedReader(new InputStreamReader(ps.getInputStream()));
            BufferedReader bre = new BufferedReader(new InputStreamReader(ps.getErrorStream()));
            String line="";
            String out="";
            while ((line=bri.readLine())!=null) {
               out+=line+"\n";
            }
            while ((line = bre.readLine()) != null) {
               out+=line+"\n";
            }
            bre.close();
            ps.waitFor();
            if (!out.toString().trim().isEmpty()) {
                  System.err.println("RNAfold says: "+out);
              }
        }
        catch(Exception e){
            e.printStackTrace();
        }
    }    
    
    /**
     * Return only one line
     * @param cmd
     * @return 
     */
    public static String executeLinuxCommandGetLine (String cmd){

        try {
            ProcessBuilder pb = new ProcessBuilder("bash", "-c", cmd);
            pb.redirectErrorStream(true); // use this to capture messages sent to stderr
            Process shell = pb.start();
            InputStream is = shell.getInputStream(); // this captures the output from the command
            int shellExitStatus = shell.waitFor(); 
            
            StringBuilder response = new StringBuilder();
              int value = 0;
              boolean active = true;
              while (active) {
                    value = is.read();
                    if (value == -1) {
                        throw new IOException("End of Stream");
                    } else if (value != '\n') {
                        response.append((char) value);
                        continue;
                    } else {
                        active = false;
                    }
                }
                return response.toString();
            
        } catch (Exception e) {
            return "";
        }
    }
    
     /**
     * Execute windows command with cmd
     * @param cmd
     */
    public static String getcmdpath() {
        Runtime runtime = Runtime.getRuntime();
        String output="";
        try{
            //Process ps=pb.start();
            Process ps=runtime.exec("cmd /c echo %cd%");
            
            BufferedReader bri=new BufferedReader(new InputStreamReader(ps.getInputStream()));
            BufferedReader bre = new BufferedReader(new InputStreamReader(ps.getErrorStream()));
            String line="";
            while ((line=bri.readLine())!=null) {
               //System.out.println(line);
               output=output+line+"\n";
            }
            while ((line = bre.readLine()) != null) {
               //System.out.println(line);
            }
            bre.close();
            ps.waitFor();
        }
        catch(Exception e){
            e.printStackTrace();
        }
        return output;
    }    

    
    /**
     * print a message in red
     * @param s 
     */
    public static void error(String s){
        System.err.println(s);
        System.exit(-1);
    }
    
    /**
     * Complement a string of IUPAC DNA nucleotides (output A C T G only).
     *
     * @param s The string of one-letter upper or lower case IUPAC nucleotides to complement.
     * @return A complemented string, ie A->T, T->A, C->G, G->C, U->A. Case is preserved.
     */
    public static String dnaComplement( String s )
    {
            char cgene[]= s.toCharArray();	

            //complement
            int i;
            for( i=0; i < cgene.length; i ++ )
            {
                    switch( cgene[ i ])
                    {
                            case 'A': cgene[ i ]= 'T'; break;
                            case 'T': cgene[ i ]= 'A'; break;
                            case 'U': cgene[ i ]= 'A'; break;
                            case 'C': cgene[ i ]= 'G'; break;
                            case 'G': cgene[ i ]= 'C'; break;
                            case 'a': cgene[ i ]= 't'; break;
                            case 't': cgene[ i ]= 'a'; break;
                            case 'u': cgene[ i ]= 'a'; break;
                            case 'c': cgene[ i ]= 'g'; break;
                            case 'g': cgene[ i ]= 'c'; break;
                    }

            }

            return new String( cgene );	
    }
    
    /**
     * Reverse and complement a string of IUPAC DNA nucleotides (A C T G only).
     *
     * @param s The string of one-letter upper or lower case IUPAC nucleotides to reverse and complement.
     * @return The reverse-complemented string, ie A->T, T->A, C->G, G->C. Case is preserved.
     */
    public static String dnaReverseComplement( String s )
    {
            StringBuffer r= new StringBuffer( dnaComplement( s ) );
            return r.reverse().toString();
    }

        /**
     * Reverse a string of IUPAC DNA/RNA nucleotides (A C T U G only).
     *
     * @param s The string of one-letter upper or lower case IUPAC nucleotides to reverse.
     * @return The reverse string. Case is preserved.
     */
    public static String Reverse( String s )
    {
            StringBuffer r= new StringBuffer(s);
            return r.reverse().toString();
    }

    /**
     * Delete a directory
     * @param path
     * @return 
     */
    static public boolean deleteDirectory(File path) {
        if( path.exists() ) {
          File[] files = path.listFiles();
          for(int i=0; i<files.length; i++) {
             if(files[i].isDirectory()) {
               deleteDirectory(files[i]);
             }
             else {
               files[i].delete();
             }
          }
        }
       return( path.delete() );
  }

    /**
     * count lines in a file
     * @param filename
     * @return 
     */
    public static int countLines(File filename) {
        int count = 0;
        try {
            InputStream is = new BufferedInputStream(new FileInputStream(filename));
            byte[] c = new byte[1024];

            int readChars = 0;
            while ((readChars = is.read(c)) != -1) {
                for (int i = 0; i < readChars; ++i) {
                    if (c[i] == '\n') {
                        ++count;

                    }
                }
            }
        } catch (IOException iOException) {
        }
        return count;
    }
            
}
