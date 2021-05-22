
import random

# generate two random chromosomes for testing test_genome

seq = {0: 'A', 1: 'C', 2: 'G', 3: 'T'}

random.seed(1234)

def random_chrom(name):
    with open(name, 'w') as oh:
        oh.write('>%s\n' % name.replace('.fa', ''))
        lb = 0
        for _ in range(3000):
            lb += 1
            if lb >= 80:
                oh.write('\n')
                lb = 0
            b = random.randint(0, 3)
            oh.write(seq[b])
        oh.write('\n')
    
random_chrom('chr1.fa')
random_chrom('chr2.fa')
random_chrom('chrA.fa') # named chroms
