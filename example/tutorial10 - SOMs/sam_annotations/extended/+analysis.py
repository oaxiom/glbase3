
from .gene_layer_name import gene_layer_name

for k in sorted(gene_layer_name): # Make sure tests pass in gene_layer_name
    print(len(gene_layer_name[k]), k)
print()

from .gene_layer_name import gene_layer_name
from .ct_t_type import derivation_type
print("\nDerivation Types:")
res = {}
for i in derivation_type:
    if derivation_type[i] not in res:
        res[derivation_type[i]] = 0
    res[derivation_type[i]] += 1
for k in sorted(res):
    print(res[k], k)
print()

print("Derivation type, per germ layer:")
res = {}
# counts per layer:
for layer in gene_layer_name:
    res[layer] = {}
    for ct in gene_layer_name[layer]:
        if derivation_type[ct] not in res[layer]:
            res[layer][derivation_type[ct]] = 0
        res[layer][derivation_type[ct]] += 1
        
for l in res:
    for k in sorted(res[l]):
        print(l, k, res[l][k])
print() 

from .sample_to_gse import sample_to_gse
print("Studies in germ layer:")
res = {}
for dom in gene_layer_name:
    if dom not in res:
        res[dom] = []
    
    for sample in gene_layer_name[dom]:
        studs = sample_to_gse[sample]
        res[dom] += studs
        
for k in sorted(res):
    print(k, len(set(res[k])))
print()    



from .tag_counts import tag_counts

res = []
for k in tag_counts:
    res += tag_counts[k]
    
print('Total number of tags in study', sum(res))