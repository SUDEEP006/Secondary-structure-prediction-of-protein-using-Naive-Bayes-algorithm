#python PSSP7.py -input protein_data.txt

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-input', dest='infile', required=True, help='Input file of script.')
args = parser.parse_args()
infile=args.infile



def preprocess_data():
    data = open(infile,"r")   
    #Opening the text file containing sequences and secondary structures.
    
    p = data.read()
    p = list(p)
    i = 0
    while i < len(p):
        if((ord(p[i])>=ord('0') and ord(p[i])<=ord('9')) or p[i]=="'" or p[i]=="\n" or p[i]==' '):
            p.pop(i)
            continue
        if(p[i] == "?" or p[i] == "_"):
            p[i] = "U"
        if(p[i]=='S' and p[i+1]=='t' and p[i+2]=='r'):
            p[i+1]='e'
            p[i+2]='q'
            i+=2
        i+=1
    p = "".join(p)            #Joining all the strings.
    p = p.split("Seq{}=")     #Making the string a list by splitting it in the part containing the string "Seq{}=".
    p.pop(0)
    seq = []
    structure = []
    for i in range(len(p)):
        if(i%2==0):
            seq.append(p[i])   #Making a list of sequences.
        else:
            structure.append(p[i])   #Making a list of the corresponding secondary structure classes.
    return(seq, structure)

#preprocess_data()

def update_sequence(window_size):
    l = window_size
    seq, structure = preprocess_data()
    pad_len = (l-1)//2    #Length of the dummy to be added at the beginning and at the end of sequence.
    updated_Seq = seq
    for i in range(len(updated_Seq)):
        updated_Seq[i] = "X"*pad_len + updated_Seq[i] + "X"*pad_len
    return(updated_Seq)

#update_sequence(5)

def train_the_data(window_size):
    seq, structure = preprocess_data()
    updated_Seq = update_sequence(window_size)   #Updated Sequence after adding dummy.
    input1 = []    
    output1 = []
    k = 0
    for i in updated_Seq:
        for j in range(len(structure[k])):
            input1.append(i[j:(j+window_size)])   #List of Sub-Sequences of given window size.
            output1.append(structure[k][j])       #List of corresponding secondary structures of sub-sequences.
        k+=1
        
    
    limit1 = int(len(input1)*0.7)                #The index for taking 70% of the data to make it training set.
    tr_input_data = input1[0:limit1]             #list of training input data.
    test_input_data = input1[limit1:]             #list of testing input data.
    tr_output_data = output1[0:limit1]           #list of training output data.
    test_output_data = output1[limit1:]           #list of testing output data.
    
        
    return (tr_input_data, tr_output_data, test_input_data, test_output_data)

window_size = 5
train_test_data = train_the_data(window_size)

tr_input_data = train_test_data[0]
tr_output_data = train_test_data[1]
test_input_data =train_test_data[2]
test_output_data =train_test_data[3]
j = 0
i = 0
#print("Training_Input", " Class", " Testing_Input", "  Class")
while j <len(test_output_data):
    #print("   ", tr_input_data[i], "       ", tr_output_data[i], "      ", test_input_data[j], "       ", test_output_data[j])
    j+=1
    i+=1
while i < len(tr_input_data):
    #print("   ", tr_input_data[i], "       ", tr_output_data[i])
    i+=1

def class_separation(tr_input_data,tr_output_data):
    dssp_sec_list = ['G','H','I','B','E','T','S','U']
    sec_structure_dict = {}
    for i in dssp_sec_list:
        sec_structure_dict[i] = []
    #print(sec_structure_dict)
    for i in range(len(tr_output_data)):
        sec_structure_dict[tr_output_data[i]].append(tr_input_data[i])
    #print(sec_structure_dict)
    #print(dssp_sec_list)
    #print(sec_structure_dict['I'])
    return (sec_structure_dict, dssp_sec_list)

sec_structure_dict, dssp_sec_list = class_separation(tr_input_data, tr_output_data)

def prob_class1(sec_structure_dict, dssp_sec_list):
    array_amino = "ARNDCQEGHILKMFPSTWYVX"
    aa_count = {}
    for i in dssp_sec_list:
        aa_count[i] = {}
        for j in array_amino:
            aa_count[i][j] = [0 for i in range(window_size)]
    class_p = {}
    P_tables = {}
    for i in dssp_sec_list:
        P_tables[i] = {}
        class_p[i] = len(sec_structure_dict[i])/len(tr_input_data)
   
    for i in dssp_sec_list:
        for j in array_amino:
            P_tables[i][j] = [0 for i in range(window_size)]
   
    for i in dssp_sec_list:
        for j in sec_structure_dict[i]:
            for k in range(window_size):
                aa_count[i][j[k]][k]+=1
            

    for i in dssp_sec_list:
        for j in array_amino:
            for k in range(window_size):
                if(len(sec_structure_dict[i]) > 0):
                    P_tables[i][j][k] = aa_count[i][j][k]/len(sec_structure_dict[i])
                else:
                    P_tables[i][j][k] = 0
    #print(P_tables)
    #print(class_p)
    return P_tables, class_p

P_tables, class_p = prob_class1(sec_structure_dict, dssp_sec_list)

testclass1 = [[0 for i in range(len(dssp_sec_list))] for j in range(len(test_input_data))]
for i in range(len(test_input_data)):
    for j in range(len(dssp_sec_list)):
        sum = 1
        for k in range(len(test_input_data[i])):
            sum*=P_tables[dssp_sec_list[j]][test_input_data[i][k]][k]
        sum*=class_p[dssp_sec_list[j]]
        testclass1[i][j] = sum
#print(dssp_sec_list)
#print(sum)
#print(testclass1)

count = 0
#print("Testing_Input", " DSSP_Class", " Predicted_ss_Class")
for i in range(len(test_input_data)):
    #print(test_input_data[i], "            ", test_output_data[i], "           ", dssp_sec_list[testclass1[i].index(max(testclass1[i]))])
    if(test_output_data[i]==dssp_sec_list[testclass1[i].index(max(testclass1[i]))]):
        count+=1
print("The Percentage Accuracy is: ", (count/len(test_input_data))*100)