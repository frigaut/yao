include,"yao.i";

parpath="./:"+pathform(_(Y_USER,Y_SITES,Y_SITE));
tmp = find_in_path("sh6x6.par",takefirst=1,path=parpath);
if (tmp==[]) tmp=find_in_path("data/sh6x6.par",takefirst=1,path=parpath);
if (tmp==[]) tmp=find_in_path("share/yao/examples/sh6x6.par",takefirst=1,path=parpath);
if (tmp!=[]) yaopardir = dirname(tmp);
if (noneof(yaopardir)) error,"Can not find parfile example directory";
cd,yaopardir;
include,"test-all.i";
