#!/bin/python
# coding: utf-8

from optparse import OptionParser
import os
import sys
import re

'''
    Title:       LDSC.py
    Description:
    Author:       Zunpeng Liu
    Email:        zunpengATmit.edu
    Date:         11/4/2022
    Version:      0.1
    Python:       based on python 3.7
    
'''


MY_USAGE='''
    python LDSC.py -i /path/to/input -p prefix -o /path/to/output -q /path/to/qsub -c /path/to/configure'''

parser=OptionParser(MY_USAGE,version='Version 1.0')
parser.add_option('-i','--in_dir',   dest='in_dir',   type='string', help='You should sepecify the input file path, wich stored all raw data.')
parser.add_option('-p','--prefix',   dest='prefix',   type='string', help='You should specifty the prefix of prefixs, which will be used to label all processed results. For example, if the prefix was defined as "M91_HHVK7CCXY_L1", the prefix of file names will be precessed as "M91".')
parser.add_option('-o','--out_dir',  dest='out_dir',  type='string', help='You should specifty the output directory, which will store all processed results.')
parser.add_option('-q','--qsub_dir', dest='qsub_dir', type='string', help='You should specifty the qub directory, which will store all shell script derived from this pipeline.')
parser.add_option('-c','--configure_file',dest='configure_file',type='string', default='/pfs1/liuguanghui/liuzunpeng/01_Script/03_ATAC-seq/07_ATAC_seq_pipeline/config.txt',help='The configure lib file should include the absolute path of softwares used in this pipeline')
(options,args)=parser.parse_args()


''' Predefine variables '''
in_dir  = str(options.in_dir)
prefix  = str(options.prefix)
out_dir = str(options.out_dir)
qsu_dir = str(options.qsub_dir)
hg19_bed  = in_dir+"/"+prefix+".bed"

''' ensure all directory has been created, if not exist,then creat the new ones. '''
def ensure_dir(directory):
    #directory = os.path.dirname(f)
    if not os.path.exists(directory):
        os.makedirs(directory)

''' make new directions '''
ensure_dir(out_dir)
ensure_dir(qsu_dir)


''' Configure all softwares and library '''
sft={}
config= options.configure_file
with open(config,'r') as config:
    for line in config:
        #        if len(line)!=0 and not line.startswith('#'):
        if len(re.findall('=',line))==1:
            sft[re.findall(r'(\w+?)=.*',line)[0]]=re.findall(r'\w+=(.*)',line)[0]
globals().update(sft)



''' 2. Creating an annot file '''
ensure_dir(str(out_dir+'/02_ldsc'))
ensure_dir(str(qsu_dir+'/02_ldsc_annot'))

ldsc_annot_str_p1='#!/bin/sh\necho "Start creating an annot file of '+prefix+'"\n'\
+'source activate\n'\
+'conda activate ldsc\n'

ldsc_annot_str_shell_file=open(qsu_dir+'/02_ldsc_annot/ldsc_annot.'+prefix+'.sh','w')
ldsc_annot_str_shell_file.write(ldsc_annot_str_p1)
ldsc_annot_str_shell_file.close()

for chromosome in range(1,23):
    annot_file=out_dir+'/02_ldsc/'+prefix+'.chr'+str(chromosome)+'.annot.gz'
    bimfile='./zunpeng/03_Database/GWAS/LDSCORE/1000G_EUR_Phase3_plink/1000G.EUR.QC.'+str(chromosome)+'.bim'	
    iteration_str=python+' '+make_annot+' --bed-file '+hg19_bed+' --bimfile '+bimfile+' --annot-file '+annot_file  +' &\n'
    ldsc_annot_str_shell_file=open(qsu_dir+'/02_ldsc_annot/ldsc_annot.'+prefix+'.sh','a+')
    ldsc_annot_str_shell_file.write(iteration_str)
    ldsc_annot_str_shell_file.close()
    
ldsc_annot_str_p3='echo "...ing"\nwait\necho "Creating an annot file '+prefix+' has been done!" '
ldsc_annot_str_shell_file=open(qsu_dir+'/02_ldsc_annot/ldsc_annot.'+prefix+'.sh','a+')
ldsc_annot_str_shell_file.write(ldsc_annot_str_p3)
ldsc_annot_str_shell_file.close()
os.system('chmod 755 ' + qsu_dir+'/02_ldsc_annot/ldsc_annot.'+prefix+'.sh')



''' 3. Computing LD scores with an annot file.  '''
ensure_dir(str(qsu_dir+'/03_ldsc_l2'))

ldsc_l2_str_p1='#!/bin/sh\necho "Computing LD scores of '+prefix+'" && \n'\
+'source activate\n'\
+'conda activate ldsc\n'

ldsc_l2_str_shell_file=open(qsu_dir+'/03_ldsc_l2/ldsc_l2.'+prefix+'.sh','w')
ldsc_l2_str_shell_file.write(ldsc_l2_str_p1)
ldsc_l2_str_shell_file.close()

for chromosome in range(1,23):
    annot_file=out_dir+'/02_ldsc/'+prefix+'.chr'+str(chromosome)+'.annot.gz'
    out_file=out_dir+'/02_ldsc/'+prefix+'.chr'+str(chromosome)
    hapmap3_snps='./zunpeng/03_Database/GWAS/LDSCORE/hapmap3_snps/hm.'+str(chromosome)+'.snp'
    bfile='./zunpeng/03_Database/GWAS/LDSCORE/1000G_EUR_Phase3_plink/1000G.EUR.QC.'+str(chromosome)
    iteration_str=python+' '+ldsc_py+' --l2 --bfile '+bfile+' --ld-wind-cm 1.0 --annot '+annot_file+' --out '+out_file+' --thin-annot --print-snps '+hapmap3_snps+' &&\n'
    ldsc_l2_str_shell_file=open(qsu_dir+'/03_ldsc_l2/ldsc_l2.'+prefix+'.sh','a+')
    ldsc_l2_str_shell_file.write(iteration_str)
    ldsc_l2_str_shell_file.close()

ldsc_l2_str_p3='echo "...ing"\necho "Computing LD scores '+prefix+' has been done!" '
ldsc_l2_str_shell_file=open(qsu_dir+'/03_ldsc_l2/ldsc_l2.'+prefix+'.sh','a+')
ldsc_l2_str_shell_file.write(ldsc_l2_str_p3)
ldsc_l2_str_shell_file.close()
os.system('chmod 755 ' + qsu_dir+'/03_ldsc_l2/ldsc_l2.'+prefix+'.sh')

######################################
''' Src '''
def getbed(bed):
    global beds
    with open(bed,'r') as b:
        for lines in b:
            lines=lines.strip()
            if len(lines) != 0 and not lines.startswith('#'):
                if len(re.findall('=', lines)) == 1:
                    beds[re.findall(r'(\w+?)=.*', lines)[0]] = re.findall(r'\w+=(.*)', lines)[0]


beds = {}
getbed("sumstat.txt")

''' 4. calculate the heritability enrichment. '''
ensure_dir(str(out_dir+'/04_Enrichment'))
ensure_dir(str(out_dir+'/04_Enrichment/'+prefix))
ensure_dir(str(qsu_dir+'/04_Enrichment'))

l2_file=out_dir+'/02_ldsc/'+prefix+'.chr'

baselineLD='./zunpeng/03_Database/GWAS/LDSCORE/1000G_Phase3_baselineLD_v2.2_hm3_snp/baselineLD.'
weights='./zunpeng/03_Database/GWAS/LDSCORE/weights_hm3_no_hla/weights.'
plink_files='./zunpeng/03_Database/GWAS/LDSCORE/1000G_Phase3_frq/1000G.EUR.QC.'

Enrichment_str='#!/bin/sh\necho "Start calculate the heritability enrichment of '+prefix+'" && \n'\
    +'source activate\n'\
    +'conda activate ldsc\n'
    
Enrichment_str_shell_file=open(qsu_dir+'/04_Enrichment/'+prefix+'.sh','w')
Enrichment_str_shell_file.write(Enrichment_str)
Enrichment_str_shell_file.close()

for sumstats_name, sumstats in beds.items():
    result_file=out_dir+'/04_Enrichment/'+prefix+'/'+prefix+'.'+sumstats_name
    Enrichment_str=python+' '+ldsc_py+' --h2 '+sumstats+' --ref-ld-chr '+l2_file+','+baselineLD+' --w-ld-chr '+weights+' --overlap-annot --frqfile-chr '+plink_files+' --out '+result_file+' --print-coefficients\n'
    Enrichment_str_shell_file=open(qsu_dir+'/04_Enrichment/'+prefix+'.sh','a')
    Enrichment_str_shell_file.write(Enrichment_str)
    Enrichment_str_shell_file.close()
    os.system('chmod 755 ' + qsu_dir+'/04_Enrichment/'+prefix+'.sh')




