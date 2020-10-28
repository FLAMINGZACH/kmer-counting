
#Q1处理一下fasta文件
def readTextFile():
    # filename = input("输入文件名:").strip()
    filename = input("输入文件名:")
    with open(filename)as file:
        next(file)
        list = []
        for line in file:
            line = line.replace("\n", "").replace("\r", "").upper().strip()
            list.append(line)
    str = ''.join(list)
    return str
str1 = readTextFile()

#Q2计算fasta文件碱基数目
def LENGTH(seq):
    Total=int(len(seq))
    count_A=str(seq.count('A'))
    count_T=str(seq.count('T'))
    count_C=str(seq.count('C'))
    count_G=str(seq.count('G'))
    print("DNA bases total number is: "+str(Total)+"\n","A number is: "+str(count_A)+"\n","T number is: "+str(count_T)+"\n","C number is: "+str(count_C)+"\n","G number is: "+str(count_G)+"\n")
#Q3 DNA互补链
def COMPLEMENT(sequence):
    sequence = sequence.replace('A', '{A}').replace('T', '{T}').replace('C', '{C}').replace('G', '{G}')
    string = sequence.format(A='T', T='A', C='G', G='C')
    return string
#Q4 DNA-RNA
def DNA_translate_RNA(DNAseq):
    RNAseq=DNAseq.replace("T","U")
    return RNAseq
# RNA-protein
def TRANSLATE(RNAseq):
    start_codon = 'ATG'
    stop_codon = ('TAA', 'TAG', 'TGA')
    codontable ={'AUA': 'I', 'AUC': 'I', 'AUU': 'I', 'AUG': 'M',
        'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACU': 'T',
        'AAC': 'N', 'AAU': 'N', 'AAA': 'K', 'AAG': 'K',
        'AGC': 'S', 'AGU': 'S', 'AGA': 'R', 'AGG': 'R',
        'CUA': 'L', 'CUC': 'L', 'CUG': 'L', 'CUU': 'L',
        'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCU': 'P',
        'CAC': 'H', 'CAU': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGU': 'R',
        'GUA': 'V', 'GUC': 'V', 'GUG': 'V', 'GUU': 'V',
        'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCU': 'A',
        'GAC': 'D', 'GAU': 'D', 'GAA': 'E', 'GAG': 'E',
        'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGU': 'G',
        'UCA': 'S', 'UCC': 'S', 'UCG': 'S', 'UCU': 'S',
        'UUC': 'F', 'UUU': 'F', 'UUA': 'L', 'UUG': 'L',
        'UAC': 'Y', 'UAU': 'Y', 'UAA': '_', 'UAG': '_',
        'UGC': 'C', 'UGU': 'C', 'UGA': '_', 'UGG': 'W'}
    start = RNAseq.find(start_codon)
    codons = []  # Create a codon list to store codons generated from coding seq.
    for sit in range(start, len(RNAseq), 3):
        if RNAseq[sit:sit + 3] in stop_codon:
            break
        if RNAseq[sit:sit + 3] in codontable.keys():
            codons.append(RNAseq[sit:sit + 3])
    proteinsequence = ''.join([codontable[codon] for codon in codons])  # Translate condons to protein seq.
    return proteinsequence
RESULT1=(DNA_translate_RNA(str1))

# Q5calculate di-nuclotide (2-mer) 最low的办法
def DINULC(seq):
    total = len(seq)
    AACount = seq.count('AA')
    print("AA",AACount,"{:.2%}".format(AACount/total))

    ATCount = seq.count('AT')
    print("AT",ATCount,"{:.2%}".format(ATCount/total))

    AGCount = seq.count("AG")
    print("AG",AGCount,"{:.2%}".format(AGCount/total))

    ACCount = seq.count("AC")
    print("AC",ACCount,"{:.2%}".format(ACCount/total))

    TACount = seq.count("TA")
    print("TA",TACount,"{:.2%}".format(TACount/total))

    TTCount = seq.count('TT')
    print("TT",TTCount,"{:.2%}".format(TTCount/total))

    TGCount = seq.count("TG")
    print("TG",TGCount,"{:.2%}".format(TGCount/total))

    TCCount = seq.count("TC")
    print("TC",TCCount,"{:.2%}".format(TCCount/total))

    GACount = seq.count("GA")
    print("GA",GACount,"{:.2%}".format(GACount/total))

    GTCount = seq.count("GT")
    print("GT",GTCount,"{:.2%}".format(GTCount/total))

    GGCount = seq.count("GG")
    print("GG",GGCount,"{:.2%}".format(GGCount/total))

    GCCount = seq.count("GC")
    print("GC",GCCount,"{:.2%}".format(GCCount/total))

    CACount = seq.count("CA")
    print("CA",CACount,"{:.2%}".format(CACount/total))

    CTCount = seq.count("CT")
    print("CT",CTCount, "{:.2%}".format(CTCount/total))

    CGCount = seq.count("CG")
    print("CG",CGCount,"{:.2%}".format(CGCount/total))

    CCCount = seq.count("CC")
    print("CC",CCCount,"{:.2%}".format(CCCount/total))

#Q6 calculate tri-nuclotide (3-mer)
def TRINULC(seq):
    total = len(seq)
    result = {}
    for i in range(total - 2):
        tmp = seq[i: i + 3]
        try:
            result[tmp] += 1
        except Exception:
            result[tmp] = 1
    for k, v in result.items():
        print(f"{k} Count is {v / total :.2%}")

#07 calculate k-mer
def mers(seq):
    a = input("请输入k-value:")
    k = int(a)
    total = len(seq)
    result = {}
    for i in range(total-k+1):
        tmp = seq[i: i + k]
        try:
            result[tmp] += 1
        except Exception:
            result[tmp] = 1
    for key, v in result.items():
        print(f"{key} Count is {v / total :.2%}")



def main():
    print(LENGTH(str1))
    print(COMPLEMENT(str1))
    print(DNA_translate_RNA(str1))
    print(TRANSLATE(RESULT1))
    print(DINULC(str1))
    print(TRINULC(str1))
    print(mers(str1))
main()
