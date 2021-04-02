import argparse

parser = argparse.ArgumentParser(description='Process some integers.')
parser.add_argument("--input")
parser.add_argument("--output")

args = parser.parse_args()

file = open(args.input, 'r')
fileW = open(args.output, 'w')
tekst =file.read()
count=0
my_hash={}

for line in tekst.split("\n")[1:-1]:
        splitted=line.split('\t')
        if splitted[2] in my_hash.keys():
                my_hash[splitted[2]] += [int(splitted[3])]

        else:
             	my_hash.update({splitted[2]:[int(splitted[3])]})
totalMax=0
totalMean=0

fileW.write("GeneName\tMean coverage\n")

for i in my_hash.keys():
        fileW.write(i + "\t" +  str(sum(my_hash[i])/len(my_hash[i]))+ "\n")

file.close()
fileW.close()
