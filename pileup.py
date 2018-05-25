#!/usr/bin/python3
import operator
import pandas as pd
import csv
from csv import reader

csv.field_size_limit(1000000000)

# This gets the variants, defined by GATK for each position and puts them in a dictionary
variantdict={}
with open("variantalternatives.tsv") as varalt:
        parser = csv.reader(varalt, skipinitialspace=True, delimiter='\t')
        for line in parser:
                key1 = str(line[0])
                key2 = str(line[1])
                key = "{}:{}".format(key1, key2)
                variantdict[key]=line[2]
varalt.close

# The function of the magical genotyper is defined here
def magicalgenotyper():
        "This is the magical genotyper, which defines the genotypes"
        qualified = [k for k, v in ualtvaldict.items() if v >= 3]
        soaltkey = [ele[0] for ele in sorted(ualtvaldict.items(), key=operator.itemgetter(1), reverse=True)]
        soaltval = [el[1] for el in sorted(ualtvaldict.items(), key=operator.itemgetter(1), reverse=True)]
        key = [e[0] for e in sorted(altvaldict.items(), key=operator.itemgetter(1), reverse=True)]
        val = [l[1] for l in sorted(altvaldict.items(), key=operator.itemgetter(1), reverse=True)]
        listofalleles = sorted([k for k, v in ualtvaldict.items() if v > 0 ])
        if len(qualified) == 0:
                if val[0] == val[1]:
                        genotype = "./."
                else:
                        if altvaldict.get(key[0]) > 12 and altvaldict.get(key[1]) == 0:
                                genotype = "{}/{}".format(key[0],key[0])
                        elif altvaldict.get(key[0]) > 1:
                                genotype = "{}/.".format(key[0])
                        else:
                                genotype = "./."
        elif len(qualified) == 1:
                if ualtvaldict.get(soaltkey[0]) > 9:
                        if len(listofalleles) > 1 and soaltval[1] == 2 and val[0] < val[1]*2:
                                genotype = "{}/.".format(soaltkey[0])
                        else:
                                genotype = "{}/{}".format(soaltkey[0], soaltkey[0])
                elif ualtvaldict.get(soaltkey[0]) in list(range(7,10)) and len(listofalleles) == 1:
                        genotype = "{}/{}".format(soaltkey[0], soaltkey[0])
                elif ualtvaldict.get(soaltkey[0]) in list(range(7,10)) and len(listofalleles) == 2 and soaltval[1] == 1:
                        genotype = "{}/{}".format(soaltkey[0], soaltkey[0])
                else:
                        if altvaldict.get(key[0]) > 12 and altvaldict.get(key[1]) == 0:
                                genotype = "{}/{}".format(key[0],key[0])
                        else:
                                genotype = "{}/.".format(key[0])
        elif len(qualified) == 2:
                if altvaldict.get(key[0]) > altvaldict.get(key[1])*12:
                        genotype = "{}/{}".format(key[0], key[0])
                elif altvaldict.get(key[0]) > altvaldict.get(key[1])*9:
                        genotype = "{}/.".format(soaltkey[0])
                else:
                        genotype = "{}/{}".format(qualified[0], qualified[1])
        elif len(qualified) == 3:
                if altvaldict.get(key[0]) > altvaldict.get(key[1])*12:
                        genotype = "{}/{}".format(key[0], key[0])
                elif altvaldict.get(key[0]) > altvaldict.get(key[1])*9:
                        genotype = "{}/.".format(soaltkey[0])
                elif  val[1] > val[2]*2:
                        genotype = "{}/{}".format(key[0], key[1])
                else:
                        genotype = "{}/.".format(soaltkey[0])
        else:
                genotype = "{}/.".format(soaltkey[0])
        print(chrompos, genotype, ualtvaldict)
        writer.writerow([chrompos, genotype, ualtvaldict, altvaldict])
        return

# iterating through each pileupfile, annotating variants and doing the genotyping
f = open("samplelist.txt")

for sample in f:
        sample = sample.rstrip()
        outputfile = open("thing/{}.out.genotype".format(sample), "w+", newline='')
        writer = csv.writer(outputfile, delimiter="\t")
        with open("thing/{}.out.bam.out.pileup".format(sample) ) as csvfile:
                for row in reader(csvfile, delimiter="\t", quoting=csv.QUOTE_NONE):
                        bases = list(row[4].translate({ord(c): None for c in "^]$<>"}))
                        positions = (row[6].split(","))
                        posbase = str(list(zip(positions, bases)))
                        uniqposbase = str(set(zip(positions, bases)))
                        uniqalta = uniqposbase.count("a") + uniqposbase.count("A")
                        uniqaltc = uniqposbase.count("c") + uniqposbase.count("C")
                        uniqaltg = uniqposbase.count("g") + uniqposbase.count("G")
                        uniqaltt = uniqposbase.count("t") + uniqposbase.count("T")
                        uniqaltdel = uniqposbase.count("*")
                        uniqref = uniqposbase.count(".") + uniqposbase.count("','")
                        ubasevaldict = {"A":uniqalta, "C":uniqaltc, "G":uniqaltg, "T":uniqaltt, "*":uniqaltdel}
                        alta = posbase.count("a") + posbase.count("A")
                        altc = posbase.count("c") + posbase.count("C")
                        altg = posbase.count("g") + posbase.count("G")
                        altt = posbase.count("t") + posbase.count("T")
                        altdel = posbase.count("*")
                        ref = posbase.count(".") + posbase.count("','")
                        basevaldict = {"A":alta, "C":altc, "G":altg, "T":altt, "*":altdel}
                        chrompos = "{}:{}".format(str(row[0]), str(row[1]))
                        alts = (variantdict.get(chrompos)).split(",")
                        alt1 = alts[0]
                        uvalalt1 = ubasevaldict.get(alt1)
                        valalt1 = basevaldict.get(alt1)
                        if len(alts) > 1:
                                alt2 = alts[1]
                                uvalalt2 = ubasevaldict.get(alt2)
                                valalt2 = basevaldict.get(alt2)
                                if len(alts) > 2:
                                        alt3 = alts[2]
                                        uvalalt3 = ubasevaldict.get(alt3)
                                        valalt3 = basevaldict.get(alt3)
                                        if len(alts) > 3:
                                                alt4 = alts[3]
                                                uvalalt4 = ubasevaldict.get(alt4)
                                                valalt4 = basevaldict.get(alt4)
                                                if len(alts) > 4:
                                                        alt5 = alts[4]
                                                        uvalalt5 = ubasevaldict.get(alt5)
                                                        valalt5 = ubasevaldict.get(alt5)
                                                        ualtvaldict = {"0":uniqref, "1":uvalalt1, "2":uvalalt2, "3":uvalalt3, "4":uvalalt4, "5":uvalalt5}
                                                        altvaldict = {"0":ref, "1":valalt1, "2":valalt2, "3":valalt3, "4":valalt4, "5":valalt5}
                                                        magicalgenotyper()
                                                else:
                                                        ualtvaldict = {"0":uniqref, "1":uvalalt1, "2":uvalalt2, "3":uvalalt3, "4":uvalalt4}
                                                        altvaldict = {"0":ref, "1":valalt1, "2":valalt2, "3":valalt3, "4":valalt4}
                                                        magicalgenotyper()
                                        else:
                                                ualtvaldict = {"0":uniqref, "1":uvalalt1, "2":uvalalt2, "3":uvalalt3}
                                                altvaldict = {"0":ref, "1":valalt1, "2":valalt2, "3":valalt3}
                                                magicalgenotyper()
                                else:
                                        ualtvaldict = {"0":uniqref, "1":uvalalt1, "2":uvalalt2}
                                        altvaldict = {"0":ref, "1":valalt1, "2":valalt2}
                                        magicalgenotyper()
                        else:
                                ualtvaldict = {"0":uniqref, "1":uvalalt1}
                                altvaldict = {"0":ref, "1":valalt1}
                                magicalgenotyper()
        csvfile.close()
        outputfile.close()

loci = pd.read_csv("loci.txt", sep="\t")
loci.columns = ["chrompos"]
f = open("samplelist.txt")
for sample in f:
        sample = sample.rstrip()
        df = pd.read_csv("thing/{}.out.genotype".format(sample), header=None, sep="\t")
        df.columns = ["chrompos", "{}".format(sample), "valdict1", "valdict2"]
        df = df.drop("valdict1", 1)
        df = df.drop("valdict2", 1)
        loci = loci.merge(df, how="left", on="chrompos")
loci = loci.fillna("./.")
loci.to_csv("out_output.tsv", sep="\t")
print("End of the program. :)")
