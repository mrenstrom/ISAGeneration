#
# gen.py [inputFile] [run/build]
#
#  from the list of files, generate commands to :
#      find consensus sequence for each group of similar seqeunces
#      run blat in consensus sequences
#      translate blat output to either best alignment or multi-alignment
#
from subprocess import call
from sys import argv
import csv
import os
import time
#
# command strings
#
slurm1 = "#!/bin/bash"
slurm2 = "#SBATCH --job-name="
slurm3 = "#SBATCH --mail-type=FAIL       # Mail events (NONE, BEGIN, END, FAIL, ALL)"
slurm4 = "#SBATCH --mail-user=mrenstrom@fredhutch.org  # Where to send mail"
slurm5 = "#SBATCH --ntasks=1                   # Run a single task"
slurm6 = "#SBATCH --cpus-per-task=1            # Number of CPU cores per task"
slurm7 = "#SBATCH --mem=32000mb                  # Total memory limit"
slurm8 = "#SBATCH --output=log/"
slurm9 = "date;hostname;pwd"

fSrc   = '"/fh/fast/kiem_h/grp/kiemlab/PROJECTS/RIS/Master Sequence Files/"'

cmdF = "~/ISAFil1/ISAFil1"

cmd1 = "~/ISAErrorCorrect/ISAErrorCorrect"

cmdb_1 = '"/fh/fast/kiem_h/grp/kiemlab/PROJECTS/BarcodeISA/RIS/blat/blat"'
cmdb_2 = '"/fh/fast/kiem_h/grp/kiemlab/PROJECTS/BarcodeISA/RIS/blat/"'
cmdb_4 = ' -out=blast8 -tileSize=11 -stepSize=5 -ooc=/fh/fast/kiem_h/grp/kiemlab/PROJECTS/BarcodeISA/RIS/blat/'

cmdp_1 = 'python ~/gen1/blast8toAln.py'

cmdg_1 = '~/buildGlobalISA/buildGlobalISA'

#-----------------------------------------------------------------------------------------------
# main start
#-----------------------------------------------------------------------------------------------
if (len(argv)) != 3:
  print('usage: genISA.py inpitFile run/build')
  exit(0)

scriptName, inputFile ,execute = argv
#
# read files to process
#
files = {}
merges=[]
subject = ""
blatlib = ""
blatooc = ""
#
# make boolean...run files or just build
#
if (execute == "run") or (execute == "Run"):
  execute = True
elif (execute == "build") or (execute == "Build"):
  execute = False
else:
  print("execute parameter(2) must be run or build")
  exit(0)
#
# read input file
#
with open(inputFile) as csv_file:
  csv_reader = csv.reader(csv_file, delimiter=',')
  for row in csv_reader:
    if row[0][0] == '#':
      pass
    elif row[0] == "SUBJECT":
        subject = row[1]
    elif row[0] == "BLAT":
        blatlib = row[1]
        blatooc = row[2]
    elif row[0] == "MERGE":
        merges.append(row)
    else:
      key = row[0]
      files[key] = row 

#
# verify
#
if len(files) == 0:
  print('could not read file ' + inputFile)
  exit(0)

if subject == "":
  print("could not set subject...ie Z09132")
  print("there should be a line in the file 'SUBJECT,ZXXXXX'")
  exit(0)

if (blatlib == "") or (blatooc == ""):
  print("could not set blat libraries")
  print("there should be a line in def file 'BLAT,xxxxx.2bit,xxx.ooc'")
  exit(0)
#
# fix destination directories based on subject
#
fxBase = '"/fh/fast/kiem_h/grp/kiemlab/PROJECTS/BarcodeISA/RIS/"' + subject + "/filtered/"
fbBase = '"/fh/fast/kiem_h/grp/kiemlab/PROJECTS/BarcodeISA/RIS/"' + subject + "/blat_files/"
ffBase = '"/fh/fast/kiem_h/grp/kiemlab/PROJECTS/BarcodeISA/RIS/"' + subject + "/final_files/"
#
# make sure batch,log and output directories exist
#
#
if not os.path.exists('log'):
  os.makedirs('log')
if not os.path.exists('batch'):
  os.makedirs('batch')
if not os.path.exists('../filtered'):
  os.makedirs('../filtered')
if not os.path.exists('../blat_files'):
  os.makedirs('../blat_files')
if not os.path.exists('../final_files'):
  os.makedirs('../final_files')
#
# convert each line into one batch file to run the job 
#
tranFiles = []
fileID  = 1
for testRow in files.values():
  test = testRow[0]
  vector = testRow[1]
  file_base    = subject + "_" + test
  job_name     = test
  file_script  = 'batch/' + file_base + '_script.sh'
  file_all     = file_base + ".fasta"
  file_con     = file_base + "_con.fasta"
  file_log     = file_base + "_f.log"
  file_b8      = file_base + ".blast8"
  file_aln     = file_base + ".aln"
  file_fasta   = file_base + ".fasta"
  file_tran    = file_base + ".tran"
  file_aln_log = file_base + "_log.txt"
  with open(file_script,'w') as sh:
    #
    # SLURM
    #
    sh.write(slurm1+'\n')
    sh.write(slurm2+job_name)
    sh.write(slurm3+'\n')
    sh.write(slurm4+'\n')
    sh.write(slurm5+'\n')
    sh.write(slurm6+'\n')
    sh.write(slurm7+'\n')
    sh.write(slurm8+file_log+'\n')
    sh.write(slurm9+'\n')
    #---------------------------------------------------------------------------
    # ISAFilterEX command
    #---------------------------------------------------------------------------
    sh.write(cmdF + " \\\n")
    sh.write(fxBase + " \\\n") 
    sh.write(str(fileID) + " " + subject + "_" + test + " " + vector)
    fileID += 1
    #
    # input files
    #
    for srcFile in testRow[2:]:
      sh.write(" \\\n" + srcFile)
    sh.write("\n")
    sh.write("\n")
    #---------------------------------------------------------------------------
    # blat command
    #---------------------------------------------------------------------------
    sh.write(cmdb_1 + ' \\\n')
    sh.write(cmdb_2 + blatlib + ' \\\n')
    sh.write(fxBase + file_all + ' \\\n')
    sh.write(cmdb_4 + blatooc + ' \\\n')
    sh.write(fbBase + file_b8 + '\n')
    sh.write("\n")
    #
    # blast8toTranFasta.py
    #
    sh.write(cmdp_1 + ' \\\n')
    sh.write(fbBase + file_b8 + ' \\\n')
    sh.write(fbBase + file_aln + '\n')
    #---------------------------------------------------------------------------
    # ISAErrorCorrecet command
    #---------------------------------------------------------------------------
    sh.write(cmd1 + " \\\n")
    sh.write(fxBase + file_all + " \\\n")
    sh.write(fbBase + file_aln + " \\\n") 
    sh.write(fbBase + file_tran + " \\\n")
    sh.write(fbBase + file_aln_log + " \n") 
    sh.write("\n")
    #
    # remember output (tran) file
    #
    tranFiles.append(file_tran)
    #
    # check for space in batch 
    #
    if execute == True:
      while True:
        a = os.popen("squeue -u mrenstrom").read() 
        print(a)
        b = len([char for char in a if char == '\n'])
        if b < 8:break
        time.sleep(15)

  if execute == True: 
    print("Start execution of " + file_script)
    #ir = os.popen("sbatch " + file_script).read()
    call(["sbatch",file_script])
#
# wait till all jobs complete
#
if execute == True:
  while True:
    a = os.popen("squeue -u mrenstrom").read() 
    print(a)
    b = len([char for char in a if char == '\n'])
    if b < 2:break
    time.sleep(15)
#
# create metafile for buildISA
#
file_meta = 'batch/meta_build.txt'
with open(file_meta,'w') as mt:
  mt.write("subject:" + subject + "\n")
  mt.write("path:"+"/fh/fast/kiem_h/grp/kiemlab/PROJECTS/BarcodeISA/RIS/" + subject + "/" + "\n")
  for tran in tranFiles:
    mt.write(tran + "\n")
#
# create sh command to run buildGlobalISA
#
file_script = 'batch/run_build_global_script.sh'
with open(file_script,'w') as sh:
    #
    # SLURM
    #
    sh.write(slurm1+'\n')
    sh.write(slurm2+'bld_g')
    sh.write(slurm3+'\n')
    sh.write(slurm4+'\n')
    sh.write(slurm5+'\n')
    sh.write(slurm6+'\n')
    sh.write(slurm7+'\n')
    sh.write(slurm8+'bld_g.log'+'\n')
    sh.write(slurm9+'\n')
    #
    #
    # buildGlobalISA command
    #
    sh.write(cmdg_1 + " " + file_meta + " \\\n")
#
# exec 
#
if execute == True:
  call(["sbatch",file_script])
#
# finished!
#
