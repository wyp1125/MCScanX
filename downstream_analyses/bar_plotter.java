import java.io.*;
import java.lang.*;
import java.util.*;
import java.awt.*;
import javax.imageio.ImageIO;
import java.awt.image.BufferedImage;

public class bar_plotter{
Hashtable gene_feat;
Hashtable chr_len;
Vector <Integer> index;
Vector <String> seg_s1;
Vector <String> seg_t1;
Vector <String> seg_s2;
Vector <String> seg_t2;

int xdim,ydim;
Vector <String> xchr;
Vector <String> ychr;
public void read_gff(String fpath)
{
try{
  FileInputStream fstream = new FileInputStream(fpath);
  DataInputStream in = new DataInputStream(fstream);
  BufferedReader br = new BufferedReader(new InputStreamReader(in));
  gene_feat=new Hashtable();
  chr_len=new Hashtable();
  String strLine;
  while ((strLine = br.readLine()) != null)   
    {
    String[] ss=strLine.split("\t");
    gene_feat.put(ss[1],ss[0]+","+ss[2]);
    if(!chr_len.containsKey(ss[0]))
    {
    chr_len.put(ss[0],ss[2]);
    }
    else
    {
    if(Integer.parseInt(ss[2])>Integer.parseInt(chr_len.get(ss[0]).toString()))
    {
    chr_len.put(ss[0],ss[2]);
    }
    }
    }
   }
catch (Exception e)
   {  
      System.err.println("Reading gff error: " + e.getMessage());
   }
}
public void read_synteny(String fpath)
{
try{
  FileInputStream fstream = new FileInputStream(fpath);
  DataInputStream in = new DataInputStream(fstream);
  BufferedReader br = new BufferedReader(new InputStreamReader(in));
  seg_s1=new Vector <String>();
  seg_t1=new Vector <String>();
  seg_s2=new Vector <String>();
  seg_t2=new Vector <String>();
  String strLine;
  String gene1="";
  String gene2="";
  int j=0;
  int strand=1;
  int i=0;
  String[] ss;
  while ((strLine = br.readLine()) != null)
  {
  if(strLine.length()>0)
  {
  if(strLine.charAt(0)=='#')
  {
  if(strLine.indexOf("Alignment")>-1)
  {
  i++;
  if(i>1)
  {
  seg_t1.add(gene1);
  if(strand==1)
  seg_t2.add(gene2);
  else
  seg_s2.add(gene2);
  }
  j=0;
  strand=1;
  if(strLine.indexOf("minus")>-1)
  strand=0;
  }
  continue;
  }
  ss=strLine.split("\t");
  gene1=ss[1];
  gene2=ss[2];
  if(j==0)
  {
  seg_s1.add(gene1);
  if(strand==1)
  seg_s2.add(gene2);
  else
  seg_t2.add(gene2);
  }
  j++;
  }
  }
  seg_t1.add(gene1);
  if(strand==1)
  seg_t2.add(gene2);
  else
  seg_s2.add(gene2);
}
catch (Exception e)
   {
      System.err.println("Reading control file error: " + e.getMessage());
   }
}

public void read_ctl(String fpath)
{
try{
  FileInputStream fstream = new FileInputStream(fpath);
  DataInputStream in = new DataInputStream(fstream);
  BufferedReader br = new BufferedReader(new InputStreamReader(in));
  String strLine;
  String[] ss,cc;
  xchr=new Vector <String>();
  ychr=new Vector <String>();
  int i;
  strLine = br.readLine();
  ss=strLine.split("\t");
  xdim=Integer.parseInt(ss[0]);
  strLine = br.readLine();
  ss=strLine.split("\t");
  ydim=Integer.parseInt(ss[0]);
  strLine = br.readLine();
  ss=strLine.split("\t");
  cc=ss[0].split(",");
  for(i=0;i<cc.length;i++)
  {
  xchr.add(cc[i]);
  }
  strLine = br.readLine();
  ss=strLine.split("\t");
  cc=ss[0].split(",");
  for(i=0;i<cc.length;i++)
  {
  ychr.add(cc[i]);
  }
}
catch (Exception e)
   {
      System.err.println("Reading synteny error: " + e.getMessage());
   }
}
public void paint (Graphics g) {
    g.setColor(Color.white);
    g.fillRect(0,0,xdim,ydim);
    g.setColor(Color.black);
    int hmargin=(int)((float)xdim*0.05);
    int vmargin=(int)((float)ydim*0.05);
    int mid_y=(int)((float)ydim/2);
    int mid_margin=vmargin;
    int rbars=xchr.size();
    int tbars=ychr.size();
    int bar_width=10;
    Double max_len=0.0;
    int i;
    int ncells=rbars;
    if(tbars>rbars)
    ncells=tbars;
    double cell_width=(double)(xdim-2*hmargin)/(double)ncells;
    for(i=0;i<rbars;i++)
    {
    if(Double.parseDouble(chr_len.get(xchr.get(i)).toString())>max_len)
    max_len=Double.parseDouble(chr_len.get(xchr.get(i)).toString());
    }
    for(i=0;i<tbars;i++)
    {
    if(Double.parseDouble(chr_len.get(ychr.get(i)).toString())>max_len)
    max_len=Double.parseDouble(chr_len.get(ychr.get(i)).toString());
    }
    double unit=(double)(mid_y-hmargin-mid_margin)/max_len;
    int temp_h,temp_x;
    Hashtable colr=new Hashtable();
    Hashtable colg=new Hashtable();
    Hashtable colb=new Hashtable();
    Hashtable xpos=new Hashtable();
    Hashtable ypos=new Hashtable();
    int cround=(rbars-1)/7+1;
    for(i=0;i<rbars;i++)
    {
    float []temp={(float)0,(float)0,(float)0};
    int change=i%7;
    float tempa=(float)0.8*(float)(i/7)/(float)cround;
////////////////////////////////////////////////////////
    switch(change){
    case 0:
    temp[0]=1-tempa;
    break;
    case 1:
    temp[1]=1-tempa;
    break;
    case 2:
    temp[2]=1-tempa;
    break;
    case 3:
    temp[0]=temp[1]=1-tempa;
    break;
    case 4:
    temp[0]=temp[2]=1-tempa;
    break;
    case 5:
    temp[1]=temp[2]=1-tempa;
    break;
    case 6:
    temp[0]=temp[1]=temp[2]=tempa;
    }
////////////////////////////////////////////////////////
//    System.out.println(tempa);
    Color syn_col=new Color(temp[0],temp[1],temp[2]);
    colr.put(xchr.get(i),temp[0]);
    colg.put(xchr.get(i),temp[1]);
    colb.put(xchr.get(i),temp[2]);
    g.setColor(syn_col);
    temp_x=(int)((double)hmargin+(double)(i)*cell_width+cell_width/3);
    xpos.put(xchr.get(i),temp_x);
    temp_h=(int)(Double.parseDouble(chr_len.get(xchr.get(i)).toString())*unit);
    g.fillRect(temp_x,(int)(mid_y-mid_margin-temp_h),10,temp_h);
    g.setColor(Color.black);
    g.drawString(xchr.get(i).toString(),temp_x-5,mid_y-5);
    } 
    for(i=0;i<tbars;i++)
    {
    temp_x=(int)((double)hmargin+(double)(i)*cell_width+cell_width/3);
    ypos.put(ychr.get(i),temp_x);
    temp_h=(int)(Double.parseDouble(chr_len.get(ychr.get(i)).toString())*unit);
    g.drawRect(temp_x,(int)(ydim-vmargin-temp_h),10,temp_h);
    g.setColor(Color.black);
    g.drawString(ychr.get(i).toString(),temp_x-5,ydim-5);
    }
    String loc_s1,loc_t1,loc_s2,loc_t2;
    String[] ss_s1,ss_t1,ss_s2,ss_t2;
    for(i=0;i<seg_s1.size();i++)
    {
    loc_s1=gene_feat.get(seg_s1.get(i)).toString();
    loc_s2=gene_feat.get(seg_s2.get(i)).toString();
    ss_s1=loc_s1.split(",");
    ss_s2=loc_s2.split(",");
    loc_t1=gene_feat.get(seg_t1.get(i)).toString();
    loc_t2=gene_feat.get(seg_t2.get(i)).toString();
    ss_t1=loc_t1.split(",");
    ss_t2=loc_t2.split(",");
    float []fcol=new float[3];
    Color fil_col;
    int fil_start, fil_len;
    if(xpos.containsKey(ss_s1[0])&&ypos.containsKey(ss_s2[0]))
    {
    fcol[0]=Float.parseFloat(colr.get(ss_s1[0]).toString());
    fcol[1]=Float.parseFloat(colg.get(ss_s1[0]).toString());
    fcol[2]=Float.parseFloat(colb.get(ss_s1[0]).toString());
    fil_col=new Color(fcol[0],fcol[1],fcol[2]);
    g.setColor(fil_col);
    fil_start=(int)((Double.parseDouble(chr_len.get(ss_s2[0]).toString())-Double.parseDouble(ss_s2[1]))*unit);
    fil_len=(int)((Double.parseDouble(ss_t2[1])-Double.parseDouble(ss_s2[1]))*unit);
    g.fillRect(Integer.parseInt(ypos.get(ss_s2[0]).toString())+1,(int)(ydim-vmargin-fil_start),8,fil_len);
    }
    if(xpos.containsKey(ss_s2[0])&&ypos.containsKey(ss_s1[0]))
    {
    fcol[0]=Float.parseFloat(colr.get(ss_s2[0]).toString());
    fcol[1]=Float.parseFloat(colg.get(ss_s2[0]).toString());
    fcol[2]=Float.parseFloat(colb.get(ss_s2[0]).toString());
    fil_col=new Color(fcol[0],fcol[1],fcol[2]);
    g.setColor(fil_col);
    fil_start=(int)((Double.parseDouble(chr_len.get(ss_s1[0]).toString())-Double.parseDouble(ss_s1[1]))*unit);
    fil_len=(int)((Double.parseDouble(ss_t1[1])-Double.parseDouble(ss_s1[1]))*unit);
    g.fillRect(Integer.parseInt(ypos.get(ss_s1[0]).toString())+1,(int)(ydim-vmargin-fil_start),8,fil_len);
    }
    } 
  }

public static void main(String args[])
{
if(args.length<8)
{
System.out.println("Usage: java bar_plotter -g gff_file -s synteny_file -c control_file -o output_PNG_file");
System.exit(1);
}
HashMap option = new HashMap();
int i;
for(i=0;i<args.length/2;i++)
{
option.put(args[2*i],args[2*i+1]);
}
if(!option.containsKey("-g")||!option.containsKey("-s")||!option.containsKey("-c")||!option.containsKey("-o"))
{
System.out.println("Usage: java bar_plotter -g gff_file -s synteny_file -c control_file -o output_PNG_file");
System.exit(1);
}
bar_plotter proc=new bar_plotter();
proc.read_gff(option.get("-g").toString());
proc.read_synteny(option.get("-s").toString());
proc.read_ctl(option.get("-c").toString());
try  {
      BufferedImage bi = new BufferedImage(proc.xdim, proc.ydim, BufferedImage.TYPE_INT_ARGB);
      Graphics2D ig2 = bi.createGraphics();
      proc.paint(ig2);
      ImageIO.write(bi, "PNG", new File(option.get("-o").toString()));
      } 
      catch (IOException ie)
      {ie.printStackTrace();}
}
}
