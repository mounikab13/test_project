
'''
##### This code takes in input as hmmscan domtblout file and gives tab separated file as output with modified column_names

#### define a function which takes input file and gives output file with desired changes

### hmmscan gives a table with header lines as below:


#                                                                                            --- full sequence --- -------------- this domain -------------   hmm coord   ali coord   env coord
# target name        accession   tlen query name                           accession   qlen   E-value  score  bias   #  of  c-Evalue  i-Evalue  score  bias  from    to  from    to  from    to  acc description of target
#------------------- ---------- -----                 -------------------- ---------- ----- --------- ------ ----- --- --- --------- --------- ------ ----- ----- ----- ----- ----- ----- ----- ---- ---------------------


###first step in the function is to extract these three lines into different variables and thereby obtaining the column names accordingly

###second line contains the column names and third line contains corresponding '---' entry for every column 

### splitting third line gives at blank spaces gives number of columns in the hmmscan output file

### first line gives information about the values correspond to which field like hmm or full sequence or alignment or envelope

###based on these three lines we can define our column names to understand easily the output file for further analysis


'''

import sys

def sys_exit(msg, err=1):
    sys.stderr.write(msg.rstrip() + "\n")
    sys.exit(err)


def convert(input_file, output_file):
    h1 = input_file.readline()
    assert h1.startswith("# ")
    h2 = input_file.readline()
    assert h2.startswith("# target name")
    h3 = input_file.readline()
    assert h3.startswith("#---")
    columns = len(h3.split())
    assert columns == 23, columns
    if columns == 23:
        names = ["target name", "accession", "tlen", "query name",
                 "accession", "qlen", "E-value", "score", "bias",
                 "#", "of", "c-Evalue", "i-Evalue", "score", "bias",
                 "from", "to", "from", "to", "from", "to",
                 "acc", "description of target"]
        assert " ".join(h2[2:-1].split()) == " ".join(names)
        names = ["target_name", "accession", "tlen", "query_name",
                 "accession", "qlen", "full_sequence_E-value",
                 "full_sequence_score", "full_sequence_bias",
                 "dom#", "ndom",
                 "c-Evalue", "i-Evalue", "score", "bias",
                 "hmmfrom", "hmmto",
                 "alifrom", "alito",
                 "envfrom", "envto",
                 "acc", "description_of_target"]
    assert len(names) == columns
    output_file.write("#%s\n" % "\t".join(names))
    if columns != 23:
        sys_exit("Expected hmmscan output!!!")

    count = 0
    for line in input_file:
        assert line[0] != "#"
        parts = line.rstrip("\n").split(None, columns - 1)
        assert len(parts) == columns, parts
        output_file.write("\t".join(parts) + "\n")
        count += 1
    return count

inp = open("sample1.dom.txt", "r")
out = open("ltrs_data", "w")

count = convert(inp, out)

print("Converted table with %i lines" % count)
 


'''

#target_name	accession	tlen	query_name	accession	qlen	full_sequence_E-value	full_sequence_score	full_sequence_bias	dom#	ndom	c-Evalue	i-Evalue	score	bias	hmmfrom	hmmto	alifrom	alito	envfrom	envto	acc	description_of_target
LTR48	DF0000536.4	787	89a1ab60-3491-4d31-a4cd-268e4aa8d309	-	1291	0.0062	9.9	0.3	1	1	2.3e-05	0.012	9.0	0.3	87	198	822	945	801	1001	0.64	LTR48 (Long Terminal Repeat) for ERV1 endogenous retrovirus
MLT1O	DF0001026.4	542	18680ac0-d8c8-4b2e-b8ed-051cbc0511bc	-	325	0.032	8.5	6.5	1	1	0.00011	0.055	7.7	6.5	286	364	17	91	10	119	0.70	MLT1O Long Terminal Repeat for ERVL-MaLR retrotransposon
MLT-int	DF0001039.4	1735	1a4c77b6-9f90-44f9-8e37-e1d50dfc1766	-	1402	5.8e-31	101.8	0.1	1	1	2.6e-32	1.3e-30	100.6	0.1	1034	1366	948	1315	907	1329	0.75	Internal region of ERVL-MaLR retrotransposon, MLT-int subfamily
MST-int	DF0001046.4	1651	1a4c77b6-9f90-44f9-8e37-e1d50dfc1766	-	1402	3.2e-30	99.4	0.1	1	1	1.5e-31	7.5e-30	98.2	0.1	1008	1324	971	1312	950	1321	0.74	Internal region of ERVL-MaLR retrotransposon, MST-int subfamily

'''



'''

after obtaining the tab seperated file we extract only few columns which we are interested in

read the file into pandas dataframe using read_csv function with tab as separator

then sort the dataframe with respect to query_name and alifrom columns which bring all the repeated queries together with alifrom values in ascending order

'''



import pandas as pd
df = pd.read_csv("ltrs_data",sep='\t')

sorted_dt = df.sort_values(['query_name','alifrom'])


'''
after sorting the queries together extract four columns named query_name,alifrom,alito and i-Evalue

make a dictionary with query_name as key and all corresponding values as list of dictionaries with column values as items.

'''
test_dt=test_dt[['query_name','i-Evalue','alifrom','alito']]

dict1=dict(test_dt.set_index('query_name').groupby(level = 0).\
    apply(lambda x : x.to_dict(orient= 'records')))



'''
after the dictionary is created we need to define a function which takes up the dictionary as input and gives output as dictionary with query_name as key and alifrom alito coordinates and minimum i-Evalue 

'''
from collections import defaultdict
def nonoverlapping_hits(dict1):

#initiate an empty dictionary initially and then after conditions are satisfied values get appended to this dictionary

    return_dict=defaultdict(list)

#create a for loop for every query in the dict.keys() 

    for query_key in dict1.keys():

#create a list called temp_query which contains values corresponding to one query   

        temp_query=dict1[query_key]

#create a temporary dictionary for every query and also taking the first values of alifrom and alito into alifrom_to_reflist

        return_dict_temp=defaultdict(list)
        alifrom_to_reflist=[(temp_query[0]['alifrom'], temp_query[0]['alito'])]

#now the analysis is to be done for every query based on the values in temp_query

#create a temporary tuple with values of one alifrom and alito at same index into alifrom_to_temp

#also initialize a comparison counter with value 0

        for idx1 in range(len(temp_query)):    
            alifrom_to_temp=(temp_query[idx1]['alifrom'], temp_query[idx1]['alito'])
            comparison_counter=0


###now we do comparison of each and every tuple with reference value in alifrom_to_reflist
            for idx2 in range(len(alifrom_to_reflist)):
                 print('reference list', alifrom_to_reflist)
                 print('comparison sets', alifrom_to_temp, 'vs.', alifrom_to_reflist[idx2])
                 print('\n')
##the comparison is done between alito value of first tuple and alifrom of next tuple

##if alifrom of second tuple is less than alito of first tuple then those two are considered to be overlapping so we add the corresponding i-Evalues of two overlapping elements to temporary dict called return_dict_temp


                if alifrom_to_reflist[idx2][1] >= alifrom_to_temp[0]:
                    return_dict_temp[alifrom_to_reflist[idx2]].append(a[idx1]['i-Evalue'])
                    comparison_counter+=1


###if the condition is not satisfied that means the queries are not overlapping and the second tuple which is compared is now added to reference list called alifrom_to_reflist

##now the comparison is done with both the values in reference list and if condition is satisfied with one of the tuples in reference list the i-Evalue is appended to temporary dictionary 

            if comparison_counter == 0:
                alifrom_to_reflist.append((a[idx1]['alifrom'], a[idx1]['alito']))
                return_dict_temp[(a[idx1]['alifrom'], a[idx1]['alito'])].append(a[idx1]['i-Evalue'])
 
#making query_name as key return_dict is created with values of return_dict_temp
   
        return_dict[query_key]=return_dict_temp
        print(return_dict_temp)

    print(return_dict)


#now the overlapping reads are filtered and corresponding i-Evalues are     
return_dict1=defaultdict(list)

for query_key, coords_iElist in return_dict.items():
    for coords, iElist in coords_iElist.items():
        return_dict1[query_key].append([coords, min(iElist)])

print(return_dict1)
