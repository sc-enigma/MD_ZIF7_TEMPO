import os
import time as tm
from datetime import datetime

def run_md(path):
    if not os.path.exists(path):
        return 
    os.chdir(path)
    os.system('gmx_mpi mdrun -s *.tpr -v')
    while True: 
        tm.sleep(1)
        if 'confout.gro' in os.listdir('.'):
            break
    tm.sleep(10)

def write_log(path):
    log_filname = '/home/sc_enigma/Projects/MD_UIO66/system/log.txt'
    lines = []
    if os.path.exists(log_filname):
        f = open(log_filname, 'r')
        lines = [line for line in f]
        f.close()
    now = datetime.now()
    current_time = now.strftime("%m/%d/%Y, %H:%M:%S")
    lines.append(f'{current_time} gmx_mpi mdrun {path}\n')
    f = open(log_filname, 'w')
    for line in lines:
        f.write(line)
    f.close()

'''
    '/home/sc_enigma/Projects/MD_UIO66/production/zif7_tempo_lpA/prod/',
    '/home/sc_enigma/Projects/MD_UIO66/production/zif7_tempo_lpB/prod/',
'''

path_list = []
tm.sleep(3000) # !!! TEMPORARY !!!
for path in path_list:
    write_log(path)
    run_md(path)
