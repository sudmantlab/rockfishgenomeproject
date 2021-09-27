#!/usr/bin/env python3

import sys
from Bio import SeqIO
import decimal
import collections

if len(sys.argv) != 2:
    print("Please call one input file")
    sys.exit()

record_dict = SeqIO.to_dict(SeqIO.parse(sys.argv[1], 'fasta'))

long_lived_list = ["Sebastes_aleutianus","Sebastes_alutus","Sebastes_aurora",
                   "Sebastes_babcocki","Sebastes_borealis","Sebastes_crameri",
                   "Sebastes_nigrocinctus","Sebastes_ruberrimus","Sebastolobus_alascanus"]

short_lived_list = ["Helicolenus_lengerichi","Sebastes_atrovirens","Sebastes_cheni",
                   "Sebastes_dalli","Sebastes_emphaeus","Sebastes_exsul","Sebastes_hopkinsi",
                   "Sebastes_inermis","Sebastes_lentiginosus","Sebastes_minor",
                   "Sebastes_rastrelliger","Sebastes_schlegelii","Sebastes_semicinctus",
                   "Sebastes_serriceps","Sebastes_steindachneri","Sebastes_ventricosus"]

#record_dict = SeqIO.to_dict(SeqIO.parse('OrthoGroup1932.FAA', 'fasta'))

all_species =  list([*record_dict])

seq_lengths1 = [len(val1) for val1 in record_dict.values()]
uniq_seq_lengths = sorted(set(seq_lengths1))
size_uniq_seq_lengths = len(uniq_seq_lengths)
if(size_uniq_seq_lengths != 1):
    sys.exit (' The sequences must be an MSA of equal length')

#curr_long_species = list(set(all_species) & set(long_lived_list))
#curr_long_species = list(set(all_species).intersection(set(long_lived_list)))

curr_long_species = list(filter(lambda x:x in long_lived_list, all_species))
long_fore_dict = dict([(x1,record_dict[x1]) for x1 in curr_long_species])
long_fore_seq_lengths = [len(val1) for val1 in long_fore_dict.values()]
long_fore_uniq_seq_lengths = len(sorted(set(long_fore_seq_lengths)))

curr_short_species = list(filter(lambda x:x in short_lived_list, all_species))
short_fore_dict = dict([(x1,record_dict[x1]) for x1 in curr_short_species])

curr_long_background = list(filter(lambda x:x not in long_lived_list, all_species))
long_back_dict = dict([(x1,record_dict[x1]) for x1 in curr_long_background])

curr_short_background = list(filter(lambda x:x not in short_lived_list, all_species))
short_back_dict = dict([(x1,record_dict[x1]) for x1 in curr_short_background])

first_species = list(record_dict.keys())[0]
align_length1 = len(record_dict[first_species])

#list_sp_types1 = ['long_fore_dict', 'long_back_dict', 'short_fore_dict', 'short_back_dict']

print("#Type","Pos","AA","Forefreq","ForeSpnum","Backfreq","BackSpnum", sep=" ")

# Foreground for Longlived
meta_base_seq_counts_LongFore = {}
for curr_base in range(0, align_length1-1):
    curr_base_exact = curr_base + 1
    list_curr_base_seq = []
    count_curr_base_seq = {}
    meta_base_seq_counts_LongFore[curr_base_exact] = {}
    for key1, value1 in long_fore_dict.items():
        list_curr_base_seq.append(long_fore_dict[key1][curr_base])
    sort_list_curr_base_seq = []
    sort_list_curr_base_seq = sorted(list_curr_base_seq)
    if "-" in list(sort_list_curr_base_seq):
        while '-' in sort_list_curr_base_seq: sort_list_curr_base_seq.remove('-')
    if not sort_list_curr_base_seq: 
        continue
    sort_list_curr_base_seq_spnum = len(sort_list_curr_base_seq)
    uniq_list_curr_base_seq = set(sort_list_curr_base_seq)
    for curr_base_all in uniq_list_curr_base_seq:
        meta_base_seq_counts_LongFore[curr_base_exact][curr_base_all] = list([round(list_curr_base_seq.count(curr_base_all) / sort_list_curr_base_seq_spnum, 3), sort_list_curr_base_seq_spnum])
    meta_base_seq_counts_LongFore[curr_base_exact] = collections.OrderedDict(sorted(meta_base_seq_counts_LongFore[curr_base_exact].items(), key=lambda v:v[1], reverse=True))
    
# Background for Longlived
meta_base_seq_counts_LongBack = {}
for curr_base in range(0, align_length1-1):
    curr_base_exact = curr_base + 1
    list_curr_base_seq = []
    count_curr_base_seq = {}
    meta_base_seq_counts_LongBack[curr_base_exact] = {}
    for key1, value1 in long_back_dict.items():
        list_curr_base_seq.append(long_back_dict[key1][curr_base])
    sort_list_curr_base_seq = []
    sort_list_curr_base_seq = sorted(list_curr_base_seq)
    if "-" in list(sort_list_curr_base_seq):
        while '-' in sort_list_curr_base_seq: sort_list_curr_base_seq.remove('-')
    if not sort_list_curr_base_seq: 
        continue    
    sort_list_curr_base_seq_spnum = len(sort_list_curr_base_seq)
    uniq_list_curr_base_seq = set(sort_list_curr_base_seq)
    for curr_base_all in uniq_list_curr_base_seq:
        meta_base_seq_counts_LongBack[curr_base_exact][curr_base_all] = list([round(list_curr_base_seq.count(curr_base_all) / sort_list_curr_base_seq_spnum, 3), sort_list_curr_base_seq_spnum])
    meta_base_seq_counts_LongBack[curr_base_exact] = collections.OrderedDict(sorted(meta_base_seq_counts_LongBack[curr_base_exact].items(), key=lambda v:v[1], reverse=True))

# Foreground for Shortlived                
meta_base_seq_counts_ShortFore = {}
for curr_base in range(0, align_length1-1):
    curr_base_exact = curr_base + 1
    list_curr_base_seq = []
    count_curr_base_seq = {}
    meta_base_seq_counts_ShortFore[curr_base_exact] = {}
    for key1, value1 in short_fore_dict.items():
        list_curr_base_seq.append(short_fore_dict[key1][curr_base])
    sort_list_curr_base_seq = []
    sort_list_curr_base_seq = sorted(list_curr_base_seq)
    if "-" in list(sort_list_curr_base_seq):
        while '-' in sort_list_curr_base_seq: sort_list_curr_base_seq.remove('-')
    if not sort_list_curr_base_seq: 
        continue
    sort_list_curr_base_seq_spnum = len(sort_list_curr_base_seq)
    uniq_list_curr_base_seq = set(sort_list_curr_base_seq)
    for curr_base_all in uniq_list_curr_base_seq:
        meta_base_seq_counts_ShortFore[curr_base_exact][curr_base_all] = list([round(list_curr_base_seq.count(curr_base_all) / sort_list_curr_base_seq_spnum, 3), sort_list_curr_base_seq_spnum])
    meta_base_seq_counts_ShortFore[curr_base_exact] = collections.OrderedDict(sorted(meta_base_seq_counts_ShortFore[curr_base_exact].items(), key=lambda v:v[1], reverse=True))

# Background for Shortlived
meta_base_seq_counts_ShortBack = {}
for curr_base in range(0, align_length1-1):
    curr_base_exact = curr_base + 1
    list_curr_base_seq = []
    count_curr_base_seq = {}
    meta_base_seq_counts_ShortBack[curr_base_exact] = {}
    for key1, value1 in short_back_dict.items():
        list_curr_base_seq.append(short_back_dict[key1][curr_base])
    sort_list_curr_base_seq = []
    sort_list_curr_base_seq = sorted(list_curr_base_seq)
    if "-" in list(sort_list_curr_base_seq):
        while '-' in sort_list_curr_base_seq: sort_list_curr_base_seq.remove('-')
    if not sort_list_curr_base_seq: 
        continue    
    sort_list_curr_base_seq_spnum = len(sort_list_curr_base_seq)
    uniq_list_curr_base_seq = set(sort_list_curr_base_seq)
    for curr_base_all in uniq_list_curr_base_seq:
        meta_base_seq_counts_ShortBack[curr_base_exact][curr_base_all] = list([round(list_curr_base_seq.count(curr_base_all) / sort_list_curr_base_seq_spnum, 3), sort_list_curr_base_seq_spnum])
    meta_base_seq_counts_ShortBack[curr_base_exact] = collections.OrderedDict(sorted(meta_base_seq_counts_ShortBack[curr_base_exact].items(), key=lambda v:v[1], reverse=True))

#######    
### Convergent in Long lived species ###
for key2, value2 in meta_base_seq_counts_LongFore.items():
    if list(value2.keys()): # Checking for emptiness, can arise due to only gaps
        current_element1 = list(value2.keys())[0] # Storing the current AA
    else: # Otherwise go to next position
        continue
    long_fore_mainbase_freq = list(value2.values())[0][0] # Storing the frequency of highest occuring AA at position
    long_fore_mainbase_spnum = list(value2.values())[0][1] # Storing the numspecies of highest occuring AA at position
    if long_fore_mainbase_freq > 0.7: # If the foreground freq is >70%
        if current_element1 in meta_base_seq_counts_LongBack[key2]: # Check for the AA in background position
            long_back_mainbase_freq = meta_base_seq_counts_LongBack[key2][current_element1][0] # Freq of background species AA in that position 
            long_back_mainbase_spnum = meta_base_seq_counts_LongBack[key2][current_element1][1]
            if long_back_mainbase_freq <= 0.3: # If background freq is <30%
                print("Long", key2,current_element1,round(long_fore_mainbase_freq,3),long_fore_mainbase_spnum,round(long_back_mainbase_freq,3),long_back_mainbase_spnum, sep=" ")
    
    
### Convergent in Short lived ones ###
for key3, value3 in meta_base_seq_counts_ShortFore.items():
    if list(value3.keys()):
        current_element1 = list(value3.keys())[0]
    else:
        continue
    short_fore_mainbase_freq = list(value3.values())[0][0]
    short_fore_mainbase_spnum = list(value3.values())[0][1]
    if short_fore_mainbase_freq > 0.7:
        if current_element1 in meta_base_seq_counts_ShortBack[key3]:
            short_back_mainbase_freq = meta_base_seq_counts_ShortBack[key3][current_element1][0]
            short_back_mainbase_spnum = meta_base_seq_counts_ShortBack[key3][current_element1][1]
            if short_back_mainbase_freq <= 0.3:
                print("Short", key3,current_element1,round(short_fore_mainbase_freq,3),short_fore_mainbase_spnum,round(short_back_mainbase_freq,3),short_back_mainbase_spnum, sep=" ")
 
### Convergent in Long but not short ###
for key4, value4 in meta_base_seq_counts_LongFore.items():
    if list(value4.keys()):
        current_element1 = list(value4.keys())[0]
    else:
        continue
    long_fore_mainbase_freq = list(value4.values())[0][0]
    long_fore_mainbase_spnum = list(value4.values())[0][1]
    if long_fore_mainbase_freq > 0.7:
        if current_element1 in meta_base_seq_counts_ShortFore[key4]:
            short_fore_mainbase_freq = meta_base_seq_counts_ShortFore[key4][current_element1][0]
            short_fore_mainbase_spnum = meta_base_seq_counts_ShortFore[key4][current_element1][1]
            if short_fore_mainbase_freq <= 0.3:
                print("Long.vs.Short", key4,current_element1,round(long_fore_mainbase_freq,3),long_fore_mainbase_spnum,round(short_fore_mainbase_freq,3),short_fore_mainbase_spnum, sep=" ")


