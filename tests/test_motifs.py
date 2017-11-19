
# Not yet an official test case

import glob

from glbase import *

all_pwms = []

for f in glob.glob('test_data/*.motif'): # Python library for wildcard file globbing
    all_pwms.append(pwms(f, format='HOMER')) # pwms is from glbase 
    
# Look at the zeroth entry in the pwm
print((all_pwms[0]))
    
# You can add multiple lists of pwms together:
meta_pwm_list = sum(all_pwms[1:], all_pwms[0])

#print meta_pwm_list

# Use glbase to remove duplicates with the same name:
meta_pwm_list = meta_pwm_list.removeDuplicates('name')

# Save the output in HOMER-like format
meta_pwm_list.saveFASTA('test_data/meta_pwm_list.motifs')