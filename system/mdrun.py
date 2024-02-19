import os
import time as tm
from datetime import datetime

def run_md(path):
    if not os.path.exists(path):
        print('PATH DOES NOT EXIST')
        return
    os.chdir(path)
    '''if '/prod/' in path:
        os.system('ln -s ../eql/confout.gro')
        os.system('gmx_mpi grompp -f prod.mdp -c *.gro -p topol.top -o prod.tpr -maxwarn 1 -v')
        tm.sleep(10)'''
    os.system('gmx_mpi mdrun -s *.tpr -v')
    while True: 
        tm.sleep(1)
        if 'confout.gro' in os.listdir('.'):
            break
    tm.sleep(10)

def write_log(path):
    log_filname = '/home/sc_enigma/Projects/MD_ZIF7/system/log.txt'
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
    print(f'{current_time} gmx_mpi mdrun {path}\n')

'''
    '/home/sc_enigma/Projects/MD_ZIF7/production/zif7_tempo_lpA/prod/',
    '/home/sc_enigma/Projects/MD_ZIF7/production/zif7_tempo_lpB/prod/',
    '/home/sc_enigma/Projects/MD_ZIF7/production/zif7_tempo_npA/prod/',
    '/home/sc_enigma/Projects/MD_ZIF7/production/zif7_tempo_npB/prod/'
'''

path_list = [\
    '/home/sc_enigma/Projects/MD_ZIF7/production/zif7_tempo_npA/prod/',
    '/home/sc_enigma/Projects/MD_ZIF7/production/zif7_tempo_npB/prod/']

for path in path_list:
    write_log(path)
    run_md(path)
