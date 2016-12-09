import os

# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------

# remember the current path but switch temporary to the files location folder
CurrentFolder  = os.getcwd()
os.chdir('OUTPUT')

fileNames = os.popen('ls -1 estimator-*' ).read().split('\n')
fileNames.pop()

d = {}
for fname in fileNames:
    sfile   =  'state'+fname[9:]
    if (os.path.exists(sfile+'a')):
        sfsize  = int(os.path.getsize(sfile)[:-1])
        safsize = int(os.path.getsize(sfile+'a')[:-1])

        if   sfsize > safsize:
             print "State file incosistency: dat = %n > data = %n" %(sfsize, safsize)
        elif sfsize < safsize:
             print "State file incosistency: dat = %n < data = %n" %(sfsize, safsize)
             sfile = sfile + 'a'

    L = fname.split('-')
    subcommand = "--X %3d --Y % 3d --Z %3d --Ax %03d --Ay 0 --APx %03d --APy 0" %tuple(map(int, L[2:7]))
    seed = int(L[8][1:])
    subcommand += " --seed %d " %seed
    subcommand += " --state OUTPUT/%s" %sfile

    if (seed in d.keys()): d[seed] += [subcommand]
    else:                  d[seed]  = [subcommand]

#switch back to folder we used to be
os.chdir(CurrentFolder)

fileName = 'submit_z-15_03.dat'
scriptFile = open(fileName,'w')
scriptFile.write('''#!/bin/bash\n\n''')

# Create the command string and output to submit file

for seed in d.keys().sort():
    for subcommand in d[seed]:
        scriptFile.write('sqsub -q NAP_8998 -o output -e output --mpp=1G -r 7d ../../main.so %s --snake --meas 10000 --signJ -1 --beta 0.22165\n' %(subcommand))
    scriptFile.write('\n')

scriptFile.close()
os.system('chmod u+x %s'%fileName)

print '\nSubmit jobs with: ./%s\n' % (fileName)
