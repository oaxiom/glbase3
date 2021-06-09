
'''

calculate and output the mapped number of tags for appending to the 

This should drag data from the cluster generated file, not just from the RSME results.
Perhaps add another column as well?

not for import. Will write a text file with the result. 

Import that

'''

from glbase import *
import numpy
from .bad_samples import bad_samples

raw_expn = expression(filename="../../pub/rsem-genes/norm_input.tsv", format={"force_tsv": True, "skiplines": 0, "ensg": 0}, expn="column[1:]")
map_stats = genelist('mapstats_rsem.tsv', format={'force_tsv': True, 'name': 0, 'skiplines': -1, 'mapped_reads': 1, 'total_reads': 2, 'map_percent': 3})

print(map_stats)

cond_names = raw_expn.getConditionNames()

# single cell rename data
single_cell_rename_map = {
    'SS Zygote E1 C1': ['SS_Zygote_E1_C1', 'SS_Zygote_E2_C1', 'SS_Zygote_E3_C1' , 'SS_Zygote_E4_C1'],
    'SS Embryo2C early E1 C2': ['SS_Embryo2C_early_E1_C2', 'SS_Embryo2C_early_E2_C1', 'SS_Embryo2C_early_E3_C2'],
    'SS Embryo2C late E6 C2': ['SS_Embryo2C_late_E6_C2', 'SS_Embryo2C_late_E7_C1', 'SS_Embryo2C_late_E7_C1', 'SS_Embryo2C_late_E7_C2', 'SS_Embryo2C_late_E8_C2', 'SS_Embryo2C_late_E9_C1', 'SS_Embryo2C_late_E9_C2'],
    'SS Embryo2C mid E5 C2': ['SS_Embryo2C_mid_E5_C2', 'SS_Embryo2C_mid_E6_C1', 'SS_Embryo2C_mid_E6_C2', 'SS_Embryo2C_mid_E7_C1', 'SS_Embryo2C_mid_E7_C2'],
    'SS Embryo4C E1 C2': ['SS_Embryo4C_E1_C2', 'SS_Embryo4C_E1_C4', 'SS_Embryo4C_E2_C1', 'SS_Embryo4C_E2_C2', 'SS_Embryo4C_E2_C3', 'SS_Embryo4C_E2_C4', 'SS_Embryo4C_E3_C1', 'SS_Embryo4C_E3_C3', 'SS_Embryo4C_E3_C4', 'SS_Embryo4C_E4_C1', 'SS_Embryo4C_E4_C2', 'SS_Embryo4C_E4_C3', 'SS_Embryo4C_E4_C4'],
    'SS Embryo8C E1 C1': ["SS_Embryo8C_E1_C1", "SS_Embryo8C_E1_C2","SS_Embryo8C_E1_C6","SS_Embryo8C_E1_C7", "SS_Embryo8C_E1_C8", "SS_Embryo8C_E2_C1", "SS_Embryo8C_E2_C2", 
        "SS_Embryo8C_E2_C3", "SS_Embryo8C_E2_C4", "SS_Embryo8C_E2_C8","SS_Embryo8C_E5_C1","SS_Embryo8C_E5_C2","SS_Embryo8C_E5_C3","SS_Embryo8C_E5_C6",
        "SS_Embryo8C_E5_C8","SS_Embryo8C_E8_C1","SS_Embryo8C_E8_C2","SS_Embryo8C_E8_C6","SS_Embryo8C_E8_C7","SS_Embryo8C_E8_C8"],
    'SS Morula E1 C10': ['SS_Morula_E1_C10', 'SS_Morula_E1_C11', 'SS_Morula_E1_C14', 'SS_Morula_E1_C15', 'SS_Morula_E1_C2', 'SS_Morula_E1_C3', 'SS_Morula_E1_C4', 'SS_Morula_E1_C5', 'SS_Morula_E1_C6', 'SS_Morula_E1_C7',
        'SS_Morula_E1_C8', 'SS_Morula_E1_C9', 'SS_Morula_E4_C1', 'SS_Morula_E4_C3', 'SS_Morula_E4_C4', 'SS_Morula_E4_C5', 'SS_Morula_E4_C6', 'SS_Morula_E4_C7', 'SS_Morula_E5_C1', 'SS_Morula_E5_C10',
        'SS_Morula_E5_C11', 'SS_Morula_E5_C12', 'SS_Morula_E5_C13', 'SS_Morula_E5_C2', 'SS_Morula_E5_C5', 'SS_Morula_E5_C6', 'SS_Morula_E5_C7', 'SS_Morula_E5_C8', 'SS_Morula_E5_C9', 'SS_Morula_E6_C1',
        'SS_Morula_E6_C12', 'SS_Morula_E6_C2', 'SS_Morula_E6_C3', 'SS_Morula_E6_C4', 'SS_Morula_E6_C6', 'SS_Morula_E6_C7', 'SS_Morula_E6_C8'],
    # blastocyst early TE
    'SS Early blastocyst E4 C8': ['SS_Early_blastocyst_E4_C8', 'SS_Early_blastocyst_E2_C5', 'SS_Early_blastocyst_E2_C4', 'SS_Early_blastocyst_E2_C3', 'SS_Early_blastocyst_E2_C2', 'SS_Early_blastocyst_E3_C13', 
        'SS_Early_blastocyst_E2_C7', 'SS_Early_blastocyst_E4_C9', 'SS_Early_blastocyst_E4_C18', 'SS_Early_blastocyst_E4_C13', 'SS_Early_blastocyst_E2_C1'],
    # blastocyst mid TE
    'SS Mid blastocyst E2 C17': ['SS_Mid_blastocyst_E2_C17', 'SS_Mid_blastocyst_E3_C7', 'SS_Mid_blastocyst_E3_C8', 'SS_Mid_blastocyst_E2_C12', 'SS_Mid_blastocyst_E3_C6', 'SS_Mid_blastocyst_E1_C14', 'SS_Mid_blastocyst_E3_C5',
        'SS_Mid_blastocyst_E3_C4', 'SS_Mid_blastocyst_E3_C2', 'SS_Mid_blastocyst_E1_C15', 'SS_Mid_blastocyst_E1_C10', 'SS_Mid_blastocyst_E3_C23', 'SS_Mid_blastocyst_E3_C18'],
    # blastocyst late TE
    'SS Late blastocyst E2 C8': ['SS_Late_blastocyst_E2_C8', 'SS_Late_blastocyst_E2_C3', 'SS_Late_blastocyst_E1_C26', 'SS_Late_blastocyst_E1_C9', 'SS_Late_blastocyst_E1_C19', 'SS_Late_blastocyst_E1_C16'],
    # blastocyst early ICM
    'SS Early blastocyst E2 C16': ['SS_Early_blastocyst_E2_C16', 'SS_Early_blastocyst_E3_C6', 'SS_Early_blastocyst_E3_C4', 'SS_Early_blastocyst_E3_C1', 'SS_Early_blastocyst_E3_C3', 'SS_Early_blastocyst_E3_C2',
        'SS_Early_blastocyst_E2_C22', 'SS_Early_blastocyst_E2_C19', 'SS_Early_blastocyst_E4_C16', 'SS_Early_blastocyst_E4_C12', 'SS_Early_blastocyst_E3_C9', 'SS_Early_blastocyst_E4_C17', 
        'SS_Early_blastocyst_E4_C14', 'SS_Early_blastocyst_E2_C17'],
    # blastocyst mid ICM
    'SS Mid blastocyst E2 C9': ['SS_Mid_blastocyst_E2_C9', 'SS_Mid_blastocyst_E3_C13', 'SS_Mid_blastocyst_E1_C20', 'SS_Mid_blastocyst_E1_C13', 'SS_Mid_blastocyst_E1_C23', 'SS_Mid_blastocyst_E2_C6', 'SS_Mid_blastocyst_E2_C4', 'SS_Mid_blastocyst_E2_C3',
        'SS_Mid_blastocyst_E2_C1', 'SS_Mid_blastocyst_E3_C11', 'SS_Mid_blastocyst_E2_C2', 'SS_Mid_blastocyst_E1_C8', 'SS_Mid_blastocyst_E1_C24', 'SS_Mid_blastocyst_E1_C11', 'SS_Mid_blastocyst_E1_C9',
        'SS_Mid_blastocyst_E1_C12', 'SS_Mid_blastocyst_E2_C15', 'SS_Mid_blastocyst_E2_C18', 'SS_Mid_blastocyst_E2_C16', 'SS_Mid_blastocyst_E2_C5', 'SS_Mid_blastocyst_E2_C14', 'SS_Mid_blastocyst_E2_C24',
        'SS_Mid_blastocyst_E2_C13', 'SS_Mid_blastocyst_E2_C23', 'SS_Mid_blastocyst_E1_C17', 'SS_Mid_blastocyst_E2_C10', 'SS_Mid_blastocyst_E1_C19', 'SS_Mid_blastocyst_E1_C5'],
    # blastocyst PrE
    'SS Late blastocyst E2 C2': ['SS_Late_blastocyst_E2_C2', 'SS_Late_blastocyst_E2_C9', 'SS_Late_blastocyst_E1_C8', 'SS_Late_blastocyst_E1_C7', 'SS_Late_blastocyst_E2_C14', 'SS_Late_blastocyst_E1_C27', 
        'SS_Late_blastocyst_E1_C10', 'SS_Late_blastocyst_E1_C24', 'SS_Late_blastocyst_E1_C11']
    }
all_keys = sum(list(single_cell_rename_map.values()), [])

biggest = ('', 0)
smallest = ('', 0)

res = {}
mapped_reads = {}
total_reads = {}
perc_mapped = {}
for c in cond_names:
    if c in bad_samples:
        continue

    if c in all_keys:
        for k in single_cell_rename_map:
            if c in single_cell_rename_map[k]:
                c_head = k
                break
    elif '_rp' in c:
        c_head = " ".join(c.split("_")[:-1])
    else:
        c_head = c.replace('_1', '').replace('_2', '')
        c_head = " ".join(c_head.split("_"))
    # Some other specific label trimming
    c_head = c_head.replace(' DFN', '').replace(' KFN-res1', '').replace(' KFN-res2', '').replace('like KFN CX', '').replace(' KFN', '')
    c_head = c_head.replace('2C embryo', 'X2C embryo').replace('4C embryo', 'X4C embryo').replace('8C embryo', 'X8C embryo') # R compatability
    c_head = c_head.replace('X8C embryo SS', 'SS Embryo8C E1 C1')
    c_head = c_head.replace('olfactorybulb', 'olfactory bulb') # Some stragglers
    if c_head == 'oocyte': c_head = 'oocytes'
    c_head = c_head.replace('embryo 2C', 'X2C embryo')
    c_head = c_head.replace('spinal cord', 'spinalcord')
    c_head = c_head.replace('ESC E14', 'ESC')
    c_head = c_head.replace('SkMuscle', 'skeletal muscle')
    c_head = c_head.replace('ESC groundstate line2', 'ESC groundstate line1')
    c_head = c_head.replace('ESC groundstate line3', 'ESC groundstate line1')
    c_head = c_head.replace('ESC groundstate line4', 'ESC groundstate line1')
    c_head = c_head.replace('CD8T cells', 'CD8T')
    c_head = c_head.replace('CD8T untre', 'CD8T')
    if c_head == 'granulocyte': c_head = 'granulocytes'
    if c_head == 'Th2': c_head = 'CD4T Th2'
    if c_head == 'medullary thymic epithelial cell': c_head = 'medullary thymic epithelial'
    if c_head == 'cortical thymic epithelial cell': c_head = 'cortical thymic epithelial'
    
    if c_head not in res:
        res[c_head] = []
        mapped_reads[c_head] = []
        total_reads[c_head] = []
        perc_mapped[c_head] = []
    num_reads = sum(raw_expn[c])
    res[c_head].append(num_reads)
    f = map_stats._findDataByKeyLazy(key='name', value="Mm_%s" % c)
    if not f:
        print('ERROR! sample %s is missing in mapstats' % c)
        continue
    mapped_reads[c_head].append(map_stats._findDataByKeyLazy(key='name', value="Mm_%s" % c)['mapped_reads'])
    total_reads[c_head].append(map_stats._findDataByKeyLazy(key='name', value="Mm_%s" % c)['total_reads'])
    perc_mapped[c_head].append(map_stats._findDataByKeyLazy(key='name', value="Mm_%s" % c)['map_percent'])

oh = open('tag_counts.py', 'w')
oh.write('\n\ntotal_tag_counts = {\n')
for c in sorted(res):
    oh.write('    "%s": [%s],\n' % (c, ",".join([str(int(i)) for i in total_reads[c]])))
oh.write('    }\n')
    
oh.write('\n# Auto-generated code by get_tag_counts.py\n\n')
oh.write('rsem_tag_counts = {\n')
for c in sorted(res):
    oh.write('    "%s": [%s],\n' % (c, ",".join([str(int(i)) for i in res[c]])))
oh.write('    }\n')

oh.write('\n\nrsem_tag_counts_average = {\n')
for c in sorted(res):
    oh.write('    "%s": %s,\n' % (c, numpy.average([int(i) for i in res[c]])))
oh.write('    }\n')

oh.write('\n\nbowtie2_tag_counts = {\n')
for c in sorted(res):
    oh.write('    "%s": [%s],\n' % (c, ",".join([str(int(i)) for i in mapped_reads[c]])))
oh.write('    }\n')

oh.write('\n\nbowtie2_tag_counts_average = {\n')
for c in sorted(res):
    oh.write('    "%s": %s,\n' % (c, numpy.average([int(i) for i in mapped_reads[c]])))
oh.write('    }\n')

oh.write('\n\nbowtie2_tag_counts_percent_map = {\n')
for c in sorted(res):
    oh.write('    "%s": [%s],\n' % (c, ",".join([str(int(i)) for i in perc_mapped[c]])))
oh.write('    }\n')

oh.write('\n\nbowtie2_tag_counts_percent_map_average = {\n')
for c in sorted(res):
    oh.write('    "%s": %s,\n' % (c, numpy.average([int(i) for i in perc_mapped[c]])))
oh.write('    }\n')
oh.close()

# Now would be a good time to write out the number of replicates per sample key (i.e. len of 
oh = open('num_replicates.py', 'w')
oh.write('\n# Auto-generated code by get_tag_counts.py\n\nnum_replicates = {\n')
for c in sorted(res):
    oh.write('    "%s": %s,\n' % (c, len(res[c])))
oh.write('    }\n')
oh.close()

# Do some testing on the contents of res to make sure it matches all the entries in 
from .ct_t_type import derivation_type

cond_names = list(derivation_type.keys()) # This is tested as correct

print('\nFound in "res" but not in derivation_type:')
for c in res:
    if c not in cond_names:
        print('! %s' % (c,))
        
print('\nFound in "derivation_type" but not in res:')
for c in cond_names:
    if c not in res:
        print('! %s' % (c,))

