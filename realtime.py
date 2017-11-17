"""

The new realtime system, basically a set of data preprocessors for expression objects

Functions take a bunch of data and spit out an glbase expression object.

"""

import csv, numpy, copy
from . import config
from .genelist import genelist
from .expression import expression

valid_methods = ["Buganim", None, "Mean"]

def process_biomark(expression_sheets, sample_descriptor_tables, method, prefix=None, omit_samples=None, control_genes=None, missing_data_ct=35, 
    convert_under_scores=False, use_pass_fail=True, use_only_prexix=False, cull_samples_with_low_control_genes=20):
    """
    **Purpose**
        Load in some BioMark data and organise it into an expression object.
        
    **Arguments**
        expression_sheets (Required)
            A list of filenames containing the expression Ct measurements, 
        
            It should look like this (a csv file):
    
            Chip Run Info,D:\PCR DATA\LQH\lqh-0826\ChipRun.bml,1131312108,48.48 (113x),GE 48x48 Standard v1,ROX,FAM-MGB,8/26/2014 3:44:20 PM,01:28:46,BIOMARKHD198
            Application Version,4.0.1
            Application Build,20130423.0820
            Export Type,Table Results
            Quality Threshold,0.65
            Baseline Correction Method,Linear (Derivative)
            Ct Threshold Method,Auto (Global)


            Experiment Information,Experiment Information,Experiment Information,Experiment Information,Experiment Information,Experiment Information,FAM-MGB,FAM-MGB,FAM-MGB,FAM-MGB,User
            Chamber,Sample,Sample,Sample,FAM-MGB,FAM-MGB,Ct,Ct,Ct,Ct,Defined
            ID,Name,Type,rConc,Name,Type,Value,Quality,Call,Threshold,Comments
            S48-A01,S48,Unknown,1,A01,Test,14.149687023483,0.977027260905451,Pass,0.00542578158304924,
            S48-A02,S48,Unknown,1,A02,Test,13.2827171910875,0.978859410170243,Pass,0.00542578158304924,
            ...
        
        sample_descriptor_tables (Required, or None)
            If set to None, then it will attempt to collect the sample names and primer names from inside the 
            expression_sheets
        
            A set of descriptors that explain what the sample and primer names are. 
            
            If the primer and sample names are correctly labelled in the expression_sheets files,
            then set this to None and it will collect that data from the expression_sheets
            
            It should look like this (a tsv file):
            
            Sample Name Cell Name
            S01 Cell_1
            S02 Cell_2
            S03 Cell_3
            ...
            FAM-MGB Name    Primer Name
            A01 MIXL1
            A02 OTX2
            A03 SNAI1
            ...
    
            Note that the expression_sheets and samples_descriptor_tables are paired, and 
            
            len(expression_sheets) == len(samples_descriptor_tables)
                    
        missing_data_ct (Optional, default=35)
            A Ct value to input into missing data values.
        
        cull_samples_with_low_control_genes (Optional, default=20)
            remove samples where the Ct in the control_genes is >20 (by default)
            
            If you have more than 1 control_gene, both genes must pass this threshold.
            
        prefix (Optional, default=None)
            add a prefix to the 'sample names' to help in discriminating samples
       
        omit_samples (Optional, default=None)
            send a list of names to exclude from the analysis. Testing is done using 'in'
            so you can also add prefix's etc. Designed so you can omit samples such as
            'negative control' etc.
        
        convert_under_scores (Optional, default=False)
            convert underscores i.e. '_' in the sample names to spaces
        
        method (Required)
            The normalisation method required.
            
            Currently three methods are implemented:
            
            'Mean' - just mean the replicates.
            'Buganim' - implements the normalisation procedure from Buganim et al., 2012; Cell.
            None - This is primarily a debug routine to see all of the values in the Ct. You probably don't want to use this.
        
        control_genes (Required if method=="Buganim")
            You must specify one ore more control primers if the method is 'Buganim'
        
        use_pass_fail (Optional, default=True)
            Use the pass/fail calls from the BioMark.
            
    """
    if sample_descriptor_tables:
        assert len(expression_sheets) == len(sample_descriptor_tables), "expression sheets and sample descriptors must be the same length"
    assert len(expression_sheets) == len(prefix), "'prefix' is not the same length as the expression sheets"
    assert method in valid_methods, "method '%s' not found" % method
    if method == "Buganim":
        assert control_genes, "You must specify control_genes when method == 'Buganim'"
    
    if not prefix:
        prefix = [None] * len(expression_sheets)
    
    data = {}
    
    if sample_descriptor_tables:
        for expn_tab, desc_tab, prefix in zip(expression_sheets, sample_descriptor_tables, prefix):
            sample_names, primer_names = load_descriptors(desc_tab, prefix, convert_under_scores)
            expn = load_expn(expn_tab, missing_data_ct) # This already filters CtCall=Fail
        
            # First I need to combine all of the expn tables
            for s in expn:
                sam_name = sample_names[s]
                if sam_name not in data:
                    data[sam_name] = {}
                for p in expn[s]:
                    prim_name = primer_names[p]
                    if prim_name not in data[sam_name]:
                        data[sam_name][prim_name] = []
                    data[sam_name][prim_name] = data[sam_name][prim_name] + expn[s][p]
            primer_names = list(primer_names.values())
    else: # descriptors are contained in the expression_sheets
        for expn_tab, prefix in zip(expression_sheets, prefix):
            #sample_names, primer_names = load_descriptors(desc_tab, prefix, convert_under_scores)
            expn = load_expn(expn_tab, missing_data_ct, prefix=prefix, collect_sample_and_primer_names=True) # This already filters CtCall=Fail
            
            # First I need to combine all of the expn tables
            primer_names = []
            sample_names = []       
            
            for sam_name in expn:
                if sam_name not in data:
                    data[sam_name] = {}
                if sam_name not in sample_names:
                    sample_names.append(sam_name)
                    
                for prim_name in expn[sam_name]:
                    if prim_name not in data[sam_name]:
                        data[sam_name][prim_name] = []
                    if prim_name not in primer_names:
                        primer_names.append(prim_name)
                    data[sam_name][prim_name] = data[sam_name][prim_name] + expn[sam_name][prim_name]
                    
    #print data
    #print sample_names
    #print primer_names
    config.log.info("process_biomark: Found %s samples" % len(data))         
    config.log.info("process_biomark: Found %s primers" % len(primer_names))  
    
    if method:
        if method == "Mean":
            for s in data:
                for p in data[s]:
                    data[s][p] = numpy.mean(data[s][p])
            sample_names = list(data.keys())
            #print data
            primer_names = set(sum([list(data[s].keys()) for s in data], []))
            #print primer_names
            
                    
        elif method == "Buganim":
            data, sample_names, primer_names = norm_buganim(data, sample_names, primer_names, control_genes, missing_data_ct, 
                cull_samples_with_low_control_genes=cull_samples_with_low_control_genes)        
        elif method == "Buganim10":
            data, sample_names, primer_names = norm_buganim(data, sample_names, primer_names, control_genes, missing_data_ct, fudge_factor=10,
                cull_samples_with_low_control_genes=cull_samples_with_low_control_genes)     
    
    # sample_names must now be a list, primer_names must now be a set
    # package the result into the format liked by expression()
    # rows are primers.
    newl = []
    sam_order = sample_names   
    sam_order.sort()
    
    # sort out omitted samples:
    newsam_order = []
    for s in sam_order:
        if omit_samples and True in [o in s for o in omit_samples]: # omit these sample names
            continue 
        newsam_order.append(s) # Keep this sample name
    sam_order = newsam_order

    for p in primer_names: # rows
        # condition data. 
        cond_data = []
        for s in sam_order:
            if s in data: # Must be true?
                pass
                
            if p in data[s]:
                Ct = data[s][p]
            else:
                Ct = missing_data_ct # What to do with missing data?
            cond_data.append(Ct)
            
        newl.append({"name": p, "conditions": cond_data})
    
    if method is None: # Send back a genelist, an expression object will mangle the output
        newgl = []
        for p in primer_names:
            i = {"name": primer_names[p]}
            for s in sam_order:
                if p in data[s]:
                    i[sample_names[s]] = data[s][p]
                else:
                    i[sample_names[s]] = missing_data_ct
            newgl.append(i)
        gl = genelist()
        gl.load_list(newgl)
        return(gl)
        
    e = expression(loadable_list=newl, cond_names=sam_order)
    e.sort_conditions()
    return(e)
        
def load_expn(filename, missing_data_ct, prefix=None, collect_sample_and_primer_names=False):
    """
    (Internal)
    **Purpose**
        Load a BioMark experimental Ct file.
    """
    table = {}
    
    with open(filename, "rU") as fh:
        ft = csv.reader(fh)
        
        for row in ft:
            if row:
                if len(row[0]) == 7 and "S" in row[0] and "A" in row[0]: 
                    # Check it passed:
                
                    # get the Ct
                    sample = row[1]
                    primer = row[4]
                    if collect_sample_and_primer_names:
                        Ct = float(row[7])
                    else:
                        Ct = float(row[6])
                    
                    if (row[8] == "Fail") or (row[10] == "Fail"): # Sometimes seen in column 10
                        Ct = float(missing_data_ct)  
                    
                    if prefix: # prefix adding is done here if required
                        sample = '%s-%s' % (prefix, sample)
                    
                    if sample not in table:
                        table[sample] = {}
                    if primer not in table[sample]:
                        table[sample][primer] = []
                    table[sample][primer].append(Ct)
    
    return(table)
        
def load_descriptors(filename, prefix, convert_under_scores):
    """
    (Internal)
    **Purpose**
        Load a descriptor file.
    """
    with open(filename, "rU") as fh:
        ft = csv.reader(fh, dialect=csv.excel_tab)
    
        sample_names = {}
        primer_names = {}
        
        for row in ft:
            if len(row[0]) < 4 and "S" in row[0]: # Sample name
                sam_name = row[1]
                if prefix:
                    sam_name = "%s_%s" % (prefix, row[1])
                if convert_under_scores:
                    sam_name = sam_name.replace("_", " ")

                sample_names[row[0]] = sam_name
                
            if len(row[0]) < 4 and "A" in row[0]:
                primer_names[row[0]] = row[1]
            
    return(sample_names, primer_names)
    
def norm_buganim(data, sample_names, primer_names, control_genes, missing_data_ct, fudge_factor=20, cull_samples_with_low_control_genes=False):
    """
    (Internal)
    Code to process the data in the manner of Buganim
    """
    # Gene filter part 1:
    for s in data:
        for p in data[s]:
            # 1d blank replicates that have >2.0 stdev.
            stdev = numpy.std(data[s][p])
            if abs(stdev) > 2.0:
                data[s][p] = missing_data_ct
            else:
                data[s][p] = numpy.mean(data[s][p])               
    
    # Sample filter:
    mark_s_for_del = [] # mark sample for deletion:
    #inv_primer_names = {v:k for k, v in primer_names.items()}
    newsample_names = []
    for s in data:
        for g in control_genes:
            if data[s][g] >= missing_data_ct:
                mark_s_for_del.append(s)
            else:
                newsample_names.append(s) # I will keep this sample
    for s in set(mark_s_for_del):
        del data[s]
    sample_names = list(set(newsample_names)) # remove the failed samples, convert sample names into a list (not set as later I sort)
    config.log.warning("process_biomark: Deleted %s samples: %s" % (len(set(mark_s_for_del)), list(set(mark_s_for_del))))

    if cull_samples_with_low_control_genes:
        mark_s_for_del = [] # mark sample for deletion:
        newsample_names = []
        for s in data:
            for ctrl_gene in control_genes:
                if data[s][ctrl_gene] > cull_samples_with_low_control_genes:
                    mark_s_for_del.append(s)
            if s not in mark_s_for_del: # For more than 1 ctrl gene;
                newsample_names.append(s) # I will keep this sample
        for s in set(mark_s_for_del):
            del data[s]
        sample_names = list(set(newsample_names)) # remove the failed samples, convert sample names into a list (not set as later I sort)
        config.log.warning("process_biomark: Deleted %s samples: %s using cull_samples_with_low_control_genes=%s" % (len(set(mark_s_for_del)), list(set(mark_s_for_del)), cull_samples_with_low_control_genes))
        #print sample_names

    # Gene filter part 2:
    primer_names = set(primer_names) # I only need the set of primer names now
    mark_p_for_del = []
    for p in set(primer_names) - set(control_genes):
        # Can sometimes have missing values:
        all_p = []
        for s in data:
            if p in data[s]:
                all_p.append(data[s][p])
        #print s, primer_names[p], [i >= missing_data_ct for i in all_p], all_p
        if False not in [i >= missing_data_ct for i in all_p]: # If all > missing_data_ct then discard
            mark_p_for_del.append(p)
            
    config.log.warning("process_biomark: Deleted %s primers: %s" % (len(set(mark_p_for_del)), list(set(mark_p_for_del))))
    for p in set(mark_p_for_del):
        primer_names.remove(p)
        for s in data:
            if p in data[s]: # Can contain missing values
                del data[s][p]     
    
    newdata = copy.deepcopy(data)
    # Convert to ACx 
    for p in primer_names:
        for s in data:
            ctrl_expn = [data[s][ctrl] for ctrl in control_genes]
            this_expn = data[s][p]
            
            #print this_expn
            
            res = []
            for ctrl in ctrl_expn:
                #for expn in this_expn:
                ACx = fudge_factor + ctrl - this_expn # They have an arbitrary fudge factor.
                res.append(ACx) 
            #print res
                    
            newdata[s][p] = numpy.mean(res)
    return(newdata, sample_names, primer_names)