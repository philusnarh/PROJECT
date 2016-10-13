#!/usr/bin/env python
#
#python  make_ant_table.py --layout KAT7_layout.txt --tname newKAT7_Ant_Table 
#
import shutil
import sys
import os 
from pyrap.tables import table, makearrcoldesc, makecoldesc, makescacoldesc,maketabdesc
import numpy as np

def create_table (ndim=None,name=None):
    names = makescacoldesc("NAME", "aa")
    off_set = makearrcoldesc("OFFSET", 1.0, 0, [3])
    statn = makescacoldesc("STATION", "aa")
    typ = makescacoldesc("TYPE", "aa")
    mnt = makescacoldesc("MOUNT", "aa")
    flg = makescacoldesc("FLAG_ROW", False)
    pos = makearrcoldesc("POSITION", 1.0, 0, [3])
    dish_diam = makescacoldesc("DISH_DIAMETER", 12.0)
    td = maketabdesc([names,off_set,statn,typ,mnt,pos,dish_diam,flg])
    t = table(name,td,nrow=ndim)
    return t

def fill_table(tname=None,xy_file=None,use=[1,2],name="NAME",position="POSITION",mount="MOUNT",
	       offset="OFFSET",station="STATION",typE="TYPE",dish_diameter="DISH_DIAMETER",flag="FLAG_ROW"):
   
    xy = np.genfromtxt(xy_file)
    N=xy.shape[0]
    tab = create_table(N,tname)
    x = use[0]-1
    y = use[1]-1
    z = y + 1
    pos = np.ndarray([N,3])
    pos[:,0],pos[:,1],pos[:,2] = xy[:,x],xy[:,y],xy[:,z]
    
    off_set = np.ndarray([N,3])
    off_set.fill(0.0)
    
    dish_diam = np.ndarray(N)
    dish_diam.fill(25.0)
    
    flg = np.ndarray([N],dtype="bool")
    flg.fill(False)
  
    names = np.ndarray([N],dtype="|S3")
    for i in range(N):
	names[i]=np.str(i)

    mnt = np.ndarray([N],dtype="|S6")
    mnt.fill("ALT-Z")
    
    statn = np.ndarray([N],dtype="|S2")
    statn.fill("P")
    
    typ = np.ndarray([N],dtype="|S15")
    typ.fill("GROUND-BASED")
    
    tab.putcol(name,names)
    tab.putcol(position,pos)
    tab.putcol(dish_diameter,dish_diam)
    tab.putcol(mount,mnt)
    tab.putcol(offset,off_set)
    tab.putcol(station,statn)
    tab.putcol(typE,typ)
    tab.putcol(flag,flg)

def make_table(layout=None,tname=None, use=[1,2]):
   if os.path.exists(tname):
     a = input( tname+" already exists. overide (1=yes, other=no)?")
     if (a==1):
       	  shutil.rmtree(tname)
     else:
        sys.exit  
   fill_table(tname,layout,use)


import sys,os
from optparse import OptionParser
if __name__=="__main__" :
        opt = OptionParser()
        opt.set_usage('%prog options ')
        opt.set_description("Creates an antenna table and/or MS given an antenna layout \
                           (provide a configuration file to make MS, --conf and the MS name in the conf. file --msname")
        opt.add_option('-l','--layout',dest='layout',default=None,help='layout file')
	opt.add_option('-u',"--use", dest="use",default="1,2",help='x and y columns of layout, default is 1,2')
        opt.add_option('-t', '--tname',dest='tname',default=None,help='Antenna table name')
        opt.add_option('-c','--conf',dest='conf',default=None,help='configuration file name')
        opt.add_option('-m', '--msname',dest='msname',default=None,help='Name of measurement set')
        opts, args = opt.parse_args(sys.argv[1:])

        layout = opts.layout
        tname = opts.tname
	use = np.fromstring(opts.use,dtype="int",sep=",")
	make_table(layout,tname,use)
	print "#### Wrote Antenna table to: "+tname
	conf = opts.conf
	msname = opts.msname
	#
	if conf is not None:		
		os.system("makems "+conf)
		os.system('rm -f *.gds *p0.vds')
		os.system('mv *_p0 %s' %msname)
		
#	 python make_ant_table.py -l kat4.txt -t kat4_11444084505 -m KAT4_11444084505.MS -c h5to.makems.cfg
#	if  not(conf==None):
#		os.system("makems %s" %conf)
# 		if os.path.exists(msname):
#			os.system("rm -r "+msname)
# 		os.system("mv "+msname.split(".MS")[0]+"_p0 "+msname)
#		from pyrap.tables import table
#		import pylab as plt
#		tab = table(msname, readonly=False)
#		uvw = tab.getcol("UVW")
#		u = uvw[:,0]
#		v = uvw[:,1]
		#uvDist = np.sqrt(u**2 + v**2)
		#plt.hist(uvDist,bins=100,log=True,color="red")
		#plt.title(msname)
		#plt.show()
		#
		#make_ant_table.fill_table(tname="ANTENNASTABLEKUNTUNE",xy_file='vlba_kuntunse.txt')
